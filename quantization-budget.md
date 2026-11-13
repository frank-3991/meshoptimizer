# 网格资产量化误差预算分配分析

本文基于 meshoptimizer 仓库（分支 `task/14`，基线提交 `f0a684f`）的实际代码，评估量化精度预算在**位置（世界空间偏移）**、**法线/切线（方向偏差）**、**纹理坐标（采样偏移）**三类渲染误差之间的分配方式。所有结论均标注代码出处（`文件:行号`），并区分「代码可直接保证的误差上界」与「必须由调用方提供资产信息才能确定的结论」。

涉及的两层代码：

- 库层 `src/`：量化原语（`quantization.cpp`、`meshoptimizer.h` 内联函数）与顶点滤波器（`vertexfilter.cpp`），只对**已归一化到固定域的数值**负责。
- 资产层 `gltf/`：gltfpack 管线（`stream.cpp`、`mesh.cpp`、`write.cpp`、`gltfpack.cpp`），负责把带世界尺度/节点变换/材质 UV 的真实资产映射到库层的固定域，并把误差预算翻译成比特数与告警。

---

## 0. 量化原语与可推导的误差上界（库层直接保证）

这些上界只依赖比特数与量化公式，与资产内容无关，是整个预算体系的「硬通货」。

| 原语 | 代码位置 | 重建公式 / 行为 | 可推导上界 |
|---|---|---|---|
| `meshopt_quantizeUnorm` | `meshoptimizer.h:1182-1190` | `q=round(v·(2^N−1))`，重建 `q/(2^N−1)`，先 clamp 到 [0,1] | 头文件注释 `meshoptimizer.h:1066`：重建误差 ≤ `1/2^(N+1)`（归一化域内） |
| `meshopt_quantizeSnorm` | `meshoptimizer.h:1192-1202` | `q=round(v·(2^(N−1)−1))`，clamp 到 [−1,1] | 头文件注释 `meshoptimizer.h:1073`：重建误差 ≤ `1/2^N`（归一化域内） |
| `meshopt_quantizeHalf` | `quantization.cpp:14-35` | IEEE-754 fp16，下溢冲零、溢出为 inf | `meshoptimizer.h:1000-1001`：幅度范围 [6e−5, 65504]，**最大相对重建误差 5e−4** |
| `meshopt_quantizeFloat` | `quantization.cpp:37-58` | 保留 fp32 二进制、截断尾数到 N 位，非规格化冲零 | 相对误差约 `2^(−N−1)`（尾数舍入），N∈[1,23]；测试见 `demo/tests.cpp:2846-2851` |
| `meshopt_computePositionExponent` | `quantization.cpp:80-127` | 共享指数网格，`iv=round(v/2^exp)` | 保证 24 位有符号 anchor 与 `max_bits` 无符号偏移**不溢出**（见 §0.1） |
| `meshopt_encodeFilterOct` / `decodeFilterOct` | `vertexfilter.cpp:1276-1322` / `78-108` | 八面体投影 + snorm，解码重归一化 | snorm 分量误差 `1/2^N`；**单位向量角度误差无闭式上界**（见 §1.2） |
| `meshopt_encodeFilterExp` | `vertexfilter.cpp:1356-1453` | 共享/独立指数 + K 位尾数 | 受 `bits` 限制的相对误差，零值不抬高共享指数（`vertexfilter.cpp:1191-1193`） |

### 0.1 共享指数位置网格保证了什么、没保证什么

`meshopt_computePositionExponent` 的契约（`meshoptimizer.h:1018-1032`）与实现（`quantization.cpp:80-127`）只保证两件事：

1. **anchor 不溢出 24 位有符号网格**：用 `maxc`（三轴最大绝对值）约束指数，且考虑 24 位有符号网格要砍掉一位以及尾数舍入临界情形（`quantization.cpp:96-104`）。
2. **相对 anchor 的偏移不溢出 `max_bits` 无符号整数**：对每轴用保守的 `floorf(min·scale)` / `ceilf(max·scale)` 算 range，超 `max_bits` 就抬高指数（`quantization.cpp:106-124`）。

`demo/tests.cpp:2939-2973` 的 `validatePositionExponent` 在两种「最近舍入」平局方向下都断言 anchor 与 range 不溢出，证明这是**代码强制保证的容纳性（containment）**，而非误差承诺。

关键区别：该函数**不接受也不知道渲染误差目标**。它的精度旋钮是调用方传入的 `min_exp`（`meshoptimizer.h:1030`：例如 `min_exp=−10` 在米制单位下限到 1mm）。也就是说，「网格步长 `2^exp` 是否满足屏幕/世界空间误差预算」完全由调用方通过 `min_exp` 决定，库只保证在该精度下整数不溢出。`demo/main.cpp:935-959` 的 DXR 示例正是逐 meshlet 求 cluster 指数再取最大值，注释明确写出 `min_exp` 应由「最大所需精度」设定。

---

## 1. 三条受限资产流程下的预算分析

### 1.1 世界空间位置偏移受限：网格尺度与非均匀实例变换如何改变位置精度

gltfpack 默认走**定点 unorm 位置**路径（`stream.cpp:481-498`），全局共享一套量化基准：

```
offset = 全局合并 AABB 的 min          (stream.cpp:149-152)
scale  = 全局 AABB 的最大边长 extent    (stream.cpp:152)
node_scale = scale / (2^bits − 1)      (stream.cpp:168，normalized 时再乘 65535)
量化: q = quantizeUnorm((v − offset) / scale, bits)   (stream.cpp:490-492)
```

**可推导的本地空间误差上界**（`stream.cpp:157` 直接写出）：

```
error_local = scale · 0.5 / 2^(bits−1)          // 每轴，半 LSB
```

默认 `pos_bits=14`（`gltfpack.cpp:1304`），故单网格在**占据整个全局 extent** 时，本地相对误差约 `0.5/2^13 ≈ 6.1e−5`。

**网格尺度的影响（代码证据）**：全局基准意味着小网格的相对精度被稀释。`stream.cpp:160-165` 显式计算每个 extent > 1e−2 的网格的 `max_rel_error = error / bounds[i].getExtent()`，超过 5% 仅打印告警并建议 `-vpf`/`-vp N`，**不拒绝、不自动加比特**。这是「库给上界、资产层只告警」的第一个证据：误差上界可由公式推导，但是否可接受取决于调用方网格的实际尺度分布。

**非均匀实例变换的影响**：量化在节点变换**之前**的本地空间进行，运行时再由节点 `translation=offset`、`scale=node_scale`（`write.cpp:1362-1390`）还原。世界空间误差为本地误差经实例矩阵的线性映射：

- 均匀缩放 `s`：世界误差 ≈ `s · error_local`，各向同性。
- 非均匀缩放 `(sx, sy, sz)`：误差变成各向异性，轴向 `≈ sk · error_local`；若再带旋转，误差椭球被旋转到世界轴（`write.cpp:1326-1329` 把 `node_scale` 逐轴乘进实例 scale，`node.cpp:161-187` 的 `decomposeTransform` 按轴长恢复 scale，证明管线保留非均匀 scale 而不强制均匀化）。
- 偏移项 `qp.offset` 在实例化路径中被实例矩阵的旋转部分作用后再加到 translation（`write.cpp:1321-1324`），蒙皮路径同理改写 inverse bind matrix（`write.cpp:1259-1268`）。这些变换在 fp32 中进行，**变换本身不增加量化网格步长，但会把本地各向同性的半 LSB 误差放大/旋转成世界空间的各向异性误差**。

结论：库与 gltfpack 保证的是**本地空间半 LSB 上界**；世界空间上界必须由调用方给出实例矩阵的最大轴向缩放 `max|sk|` 才能推出 `error_world ≤ max|sk| · error_local`。调用方若不提供变换信息，代码无法估计世界误差。

**浮点位置路径**（`-vpf`，`stream.cpp:303-310, 477-479`）改用 `meshopt_quantizeFloat`，误差转为**相对误差**（约 `2^(−bits−1)`），不插入新节点，因此不受全局 extent 稀释，适合尺度悬殊或超大坐标场景——但这是调用方依据告警做的降级选择，不是自动决策。

### 1.2 法线方向偏差受限

法线/切线有两条路径（`stream.cpp:580-631`）：

- 默认 snorm 直量化：`encodeSnorm` → `meshopt_quantizeSnorm`（`stream.cpp:350-382`），默认 `nrm_bits=8`（`gltfpack.cpp:1306`）。
- `compressmore` 时八面体滤波：`meshopt_encodeFilterOct`（`vertexfilter.cpp:1276-1322`），先投影到 2D 再 snorm。

**可推导部分**：每个 snorm 分量的归一化域误差 ≤ `1/2^N`（`meshoptimizer.h:1073`）。八面体解码端会重归一化（`vertexfilter.cpp:96-102`），测试给出经验包络：8 位 `|decoded/127 − data| < 1e−2`、12 位 `< 1e−3`（`demo/tests.cpp:904, 933`）。

**不可纯推导、需资产测量的部分**：单位向量的**角度偏差**取决于方向在八面体映射上的位置（投影在赤道折叠处非线性，`vertexfilter.cpp:1299-1300`），不存在对所有方向统一的闭式角度上界；测试中的 1e−2/1e−3 是实测包络而非证明。法线经非均匀实例变换时还要走 inverse-transpose（`mesh.cpp:56-60, 88-107`），量化误差与变换条件数耦合，同样需要资产侧的矩阵信息才能评估最终着色误差。

资产层把法线预算接到重索引容差：`reindexTolerance` 对 normal/tangent 取 `max(2e−4, 1/2^nrm_bits)`（`mesh.cpp:774-776`），使重索引合并阈值与量化 LSB 对齐，避免「量化后本应相同的法线却保留两份顶点」。

### 1.3 纹理采样偏移受限：纹理变换对采样精度的影响

UV 走 unorm 定点路径（`stream.cpp:551-578`），基准按**材质分组**而非全局：

```
qt.offset[k] = 该材质关联网格的 UV min      (stream.cpp:274-275)
qt.scale[k]  = UV max − min（每轴独立）      (stream.cpp:276-277)
q = quantizeUnorm((uv − qt.offset)/qt.scale, tex_bits)   (stream.cpp:571-572)
```

默认 `tex_bits=12`（`gltfpack.cpp:1305`）。单轴半 LSB 归一化误差 `0.5/2^(tex_bits−1)`，乘以该轴 `qt.scale` 得到 UV 空间绝对误差。

**纹理变换（`KHR_texture_transform`）的影响**在 `writeTextureInfo`（`write.cpp:183-238`）：

- 量化的 scale/offset 被**合并进**材质的纹理变换：`transform.scale *= qt.scale/(2^bits−1)`（`write.cpp:204-205`），偏移经旋转矩阵旋转后累加（`write.cpp:199-203`）。
- 因此运行时纹理采样坐标 = 资产原有纹理变换 ∘ 量化反变换。量化误差只作用于顶点 UV 本身，纹理变换在 fp32 着色器中执行，**不额外损失精度**；但原纹理变换的 `scale` 会放大 UV 绝对误差对纹素采样的影响（重复平铺次数越多，同样的 UV 误差跨越的纹素越多）。
- 资产层的告警（`stream.cpp:231-259`）用 `computeUvArea`（`stream.cpp:96-129`，按三角形平均 UV 面积）把绝对误差归一化为 `error/sqrt(uvarea)`，超过 10% 才告警（`stream.cpp:258-259`）。**这个相对误差依赖资产实际 UV 展开密度，是典型的「需资产测量」结论**；库只能给归一化域半 LSB 上界。

### 1.4 预算保证来源汇总

| 预算项 | 库代码直接保证 | 依赖调用方提供资产信息 |
|---|---|---|
| 位置本地空间半 LSB 上界 | 是（unorm/snorm 公式、`stream.cpp:157`） | — |
| 位置世界空间误差 | — | 实例矩阵最大轴向缩放、旋转（`write.cpp:1321-1329`） |
| 共享指数不溢出 anchor/`max_bits` | 是（`quantization.cpp:96-124`，测试 `tests.cpp:2939-2973`） | 精度目标 `min_exp`、cluster AABB（`main.cpp:935-959`） |
| 法线分量误差 | 是（snorm `1/2^N`） | 角度/着色误差、inverse-transpose 条件数 |
| 八面体角度包络 | 实测包络（`tests.cpp:904,933`），非闭式 | 方向分布、变换 |
| UV 归一化域误差 | 是（unorm 半 LSB） | 相对纹素误差（UV 面积、平铺次数，`stream.cpp:246-255`） |
| 预算越界处置 | 仅告警（`stream.cpp:164-165, 258-259`） | 由调用方决定加比特或转浮点 |

---

## 2. 几何共享与边界共享：独立压缩 vs 共享量化基准

### 2.1 两种共享粒度的代码事实

- **位置：文件级全局基准**。`prepareQuantizationPosition` 把所有 mesh 的 AABB 合并成一个 `Bounds`（`stream.cpp:138-153`），全部网格共用 `offset/scale/node_scale`。
- **UV：材质级基准 + 并查集**。`prepareQuantizationTexture` 用 union-find 把「同一 mesh 上用到的材质」强制归入同一量化组（`stream.cpp:186-212`），因为同一顶点缓冲只能配一套 UV 反变换；组内再合并 AABB（`stream.cpp:217-228, 263-279`）。
- **共享指数路径（库/DXR）：mesh 级共享、cluster 级取最大值**。`demo/main.cpp:939-959` 对每个 meshlet 求 `computePositionExponent` 再取 `max`，使全 mesh 共用一个 exponent，而每 cluster 只存 16 位偏移 + 24 位 anchor。

### 2.2 取舍

**独立压缩（每网格/每 cluster 独立基准）**

- 收益：基准贴合局部 AABB，小特征获得满比特相对精度；cluster 仅存偏移，存储最省（`main.cpp:975-989` 的 anchor+offset 方案）。
- 代价：相邻网格若各自量化，接缝两侧同一世界坐标被舍入到**不同网格的不同 LSB**，产生裂缝。库不掌握邻接关系，无法阻止。

**共享量化基准（gltfpack 默认）**

- 收益：接缝两侧顶点在同一 `(offset, scale, bits)` 下量化，**相同输入必然得到相同整数**（`quantizeUnorm` 是确定性纯函数，`meshoptimizer.h:1182-1190`），从机制上消除「量化本身引入」的接缝裂缝；同时多实例共享同一顶点缓冲时只存一份几何，保留实例化存储收益（`write.cpp:1302-1360` 的 `EXT_mesh_gpu_instancing`）。
- 代价：全局 extent 被最大网格撑大，小网格相对精度下降，触发 §1.1 的 5% 告警。

### 2.3 兼顾接缝连续、误差限制与存储收益的适用条件

1. **共享几何多实例 + 相邻网格共享边界**：只要各实例共享同一顶点缓冲，天然共用基准，接缝连续且几何只存一份——这是 gltfpack 实例化路径的设计前提（`gltfpack.cpp:761-777`）。
2. **网格尺度比受控**：全局基准下，最小网格 extent 与全局 extent 之比不应使 `error/extent` 越过告警阈值（5% 位置 / 10% UV）。满足时共享基准同时满足接缝与误差，且无需额外存储。
3. **尺度悬殊时分层**：对尺度差异巨大的资产，应按空间区域/材质拆分量化组（UV 已有材质级并查集范式可借鉴），或对位置改用浮点路径 `-vpf`（相对误差，不受 extent 稀释，`stream.cpp:477-479`），在组内仍保持基准共享以保接缝。
4. **cluster 级 anchor + 共享 exponent**（DXR/Nanite 式）：适用条件是 cluster 间范围能在 `max_bits`（如 16）内容纳（`quantization.cpp:117-124` 保证），用 24 位 anchor 维持全局一致性、16 位偏移拿回存储收益。`main.cpp` 注释指出 `min_exp` 需基于 mesh AABB 设定，否则不严谨（`main.cpp:936-938`）。

### 2.4 约束冲突时由哪一层拒绝或降级

代码中**没有任何层因误差超预算而硬拒绝**，处置分为三级：

- **库层**：只做 `assert` 合法性检查（如 `quantization.cpp:82-83` 的 `min_exp≥−126`、`max_bits∈[2,24]`；`vertexfilter.cpp:1279-1280` 的 bits/stride 约束），以及共享指数的溢出保护（自动抬高 exponent，`quantization.cpp:121-124`）——溢出时**降级精度而非失败**。
- **资产层（gltfpack）**：超预算只发 `stderr` 告警并建议 `-vpf`/`-vp N`/`-vtf`/`-vt N`（`stream.cpp:164-165, 258-259`），把决策权交给调用方。
- **格式/合并层**：合并可行性由结构性条件门控（`canMergeMeshes`，`mesh.cpp:218-266`：场景、节点变换、材质、skin、targets 必须一致）；带变换的节点不允许合并（`mesh.cpp:198-202`），从结构上避免把不同变换空间的几何塞进同一基准。UV 侧由并查集强制同 mesh 材质同基准（`stream.cpp:189-212`），是「为一致性主动收窄自由度」的拒绝逻辑。
- **重打包防护**：对已经量化过的输入再次打包，解析器告警「repacking may result in increased quantization error」（`parsegltf.cpp:647`），提示二次量化误差累积——同样不拒绝。

---

## 3. 极端输入、顶点合并与误差的下游传播

### 3.1 零范围坐标（degenerate AABB / 零值分量）

- **共享指数路径**：零 AABB 不抬高指数。`optlog2` 对零/非规格化返回被 `min_exp` 钳制的值（`vertexfilter.cpp:1191-1193` 注释明确「zero component should not inflate the shared exponent」）；`demo/tests.cpp:2998-2999` 验证全零 bounds 返回 `min_exp`，即零范围时保留调用方请求的最高精度而非退化为粗网格。
- **gltfpack 定点路径**：`scale==0` 时把逆 scale 直接置 0（`stream.cpp:313, 483, 502, 562-563`），量化输入变 0、输出 unorm 0，配合 `node_scale=scale/(2^bits−1)=0`（`stream.cpp:168`），除零被显式规避，退化为「全部顶点映射到 offset 一点」，不产生 NaN/Inf。
- **Exp 滤波**：`EncodeExpSeparate` 下零值继承上一个指数以保持可压缩性（`vertexfilter.cpp:1414-1416`）；`EncodeExpClamped` 把负指数钳到 0（`vertexfilter.cpp:1419-1426`），保护小值不被噪声精度污染。

### 3.2 极端数值（溢出 / Inf / NaN / 超大坐标）

- fp16：溢出 → ±inf（`quantization.cpp:28-29`），下溢 → 0（`quantization.cpp:25-26`），NaN 统一为 qNaN（`quantization.cpp:31-32`）；测试覆盖 1e20→inf、1e−20→0（`tests.cpp:2890-2895`）。
- `quantizeFloat`：保留 inf/NaN 不参与舍入，非规格化冲零（`quantization.cpp:50-54`）。
- unorm/snorm：先 clamp 再量化（`meshoptimizer.h:1186-1187, 1198-1199`），越界值被钳到端点而非回绕；正常流程中 `offset/scale` 由数据自身 AABB 导出（`stream.cpp:149-152`），数据天然落在 [0,1]，clamp 只作为防御。
- 共享指数：端点接近 24 位有符号极限时，靠 `maxc_off` 的尾数临界修正（`quantization.cpp:99-104`）与 range 边界修正（`quantization.cpp:123-124`）抬高指数；`tests.cpp:3016-3042` 覆盖 8388607/8388608、半整数 tie 等边界，证明极端幅值下以「降低一档精度」换取不溢出。
- 自定义整数属性：exp 编码仅支持到 2^23，超过则回退原始 fp32 流（`stream.cpp:741-742`），是显式的容量降级。

### 3.3 量化后顶点合并及其下游影响

**管线顺序**（`gltfpack.cpp:444-482` 与 `mesh.cpp:1324-1348`）：

```
filterStreams → generateNormals/Tangents → mergeMeshes
→ processMesh: filterBones → reindexMesh(浮点容差) → filterTriangles
             → simplify(可选) → optimizeMesh(顶点获取重排)
→ prepareQuantizationPosition/Texture → writeVertexStream(实际量化)
```

- 重索引发生在量化**之前**，容差已按量化 LSB 预设（法线 `1/2^nrm_bits`、UV `0.5/2^tex_bits`，`mesh.cpp:770-790`），用 `meshopt_generateVertexRemapCustom` 按分量差阈值合并（`mesh.cpp:814-828`）。这保证「量化后会变得相同的属性」在浮点阶段就先合并，减少量化后的冗余顶点。
- **量化本身不做几何焊接**：`writeVertexStream` 只逐顶点映射（`stream.cpp:485-495`），量化后若两个不同顶点塌缩到同一整数，管线不会再次重映射。后果分两类：
  - 对位置：塌缩可能产生退化三角形，但 `filterTriangles`（`mesh.cpp:842-865`）在量化前已运行，**量化新增的退化三角形不会被剔除**，仅占用存储/光栅化空转，不会破坏索引合法性。
  - 对 morph target：基态位置与 delta 分别量化（`stream.cpp:500-549` 的 target 分支，delta 用对称有符号编码并按 maxv 选 8/16 位），两者舍入独立，可能使「基态+delta 重建」的顶点与相邻网格接缝不再严格一致——这是共享基准也无法完全覆盖的残余风险，需资产侧测量。
- 简化器在量化前的浮点几何上运行，误差度量是二次误差（`simplifier.cpp:751-773`），`target_error` 默认 1e−2 并按网格 quality 缩放（`mesh.cpp:1344, 1396`）。量化误差不进入简化决策，因此**简化误差与量化误差是相加的两笔独立预算**，调用方需自行留余量。

### 3.4 解码后渲染链路的误差传播

- 位置：整数经节点 `scale=node_scale` + `translation=offset` 还原（`write.cpp:1372-1383`），半 LSB 本地误差按 §1.1 经实例矩阵映射为世界空间各向异性偏移，直接表现为曲面抖动/接缝。
- 法线：八面体/snorm 解码后重归一化（`vertexfilter.cpp:96-102`），方向误差经 inverse-transpose 变换影响光照；切线手性在负行列式时翻转（`mesh.cpp:114-119`），量化不改变该逻辑。
- UV：量化 UV 经合并后的纹理变换采样（`write.cpp:197-207`），误差体现为纹素级漂移，平铺/高 UV 面积下被放大。
- 权重：8 位 unorm 后用 `renormalizeWeights` 把总和补偿回 255（`stream.cpp:332-348`），注释假设每分量误差 ≤0.5、总量 ≤2，保证补偿只改最大分量——这是一个**由代码显式利用的量化误差上界**。

### 3.5 可推导上界 vs 需资产测量的结论（汇总）

**可由代码/公式直接推导（不依赖具体资产）**

- unorm 半 LSB `1/2^(N+1)`、snorm `1/2^N`、fp16 相对误差 5e−4、fp 尾数 `2^(−N−1)`（§0 表）。
- gltfpack 位置本地误差 `scale·0.5/2^(bits−1)`（`stream.cpp:157`）、UV 单轴误差 `qt.scale·0.5/2^(bits−1)`（`stream.cpp:247`）。
- 共享指数 anchor/偏移不溢出（`quantization.cpp` + `tests.cpp` 双重保证）。
- 零范围除零保护、极端值 clamp/冲零/溢出 inf 的确定性去向（§3.1-3.2）。
- 权重补偿成立的误差前提 ≤0.5/分量（`stream.cpp:339-340`）。

**必须由调用方提供资产信息或实测才能确定**

- 世界空间误差：需实例矩阵轴向缩放/旋转（§1.1）。
- 八面体法线角度误差与最终着色偏差：方向分布 + 变换条件数（§1.2）。
- UV 相对纹素误差：UV 展开面积、纹理平铺次数（`computeUvArea`，§1.3）。
- 全局基准对小网格的相对精度稀释：实际 AABB 尺度分布（告警只报不测修复，§2.4）。
- 简化误差与量化误差的叠加是否超视觉阈值：两套独立预算，需资产侧回归。
- 二次量化（重打包）累积误差：`parsegltf.cpp:647` 仅告警。

---

## 4. 可执行的预算分配方案

以下方案把上述代码事实落成一套可操作流程；默认比特数沿用 gltfpack（`gltfpack.cpp:1303-1307`：pos 14 / tex 12 / nrm 8 / col 8）。

### 步骤 1：调用方先补齐资产侧输入（库无法获取的信息）

- 每个网格的本地 AABB、所有实例矩阵（求最大轴向缩放 `Smax`、行列式/条件数）。
- 材质 UV 变换的平铺次数、纹理分辨率与允许的纹素漂移。
- 世界空间误差目标 `ε_world`（如米制 1mm）、法线角度目标、UV 纹素目标。

### 步骤 2：位置预算

1. 由世界目标反推本地目标：`ε_local = ε_world / Smax`（非均匀变换取最严轴）。
2. 共享指数/DXR 路径：设 `min_exp = floor(log2(2·ε_local))`，逐 cluster 调 `meshopt_computePositionExponent(min_exp, max_bits=16)`，全 mesh 取最大 exponent（照搬 `demo/main.cpp:939-959`）。库保证不溢出，精度由 `min_exp` 锁定。
3. gltfpack 定点路径：先按默认 14 位跑一遍，读取 `stream.cpp:160-165` 的相对误差告警；任一网格 `error/extent > 5%` 时，按 `bits ≥ ceil(log2(scale/ε_local)) + 1` 加位（`-vp N`），或对尺度悬殊资产切 `-vpf` 走相对误差浮点路径。
4. 接缝优先的资产保持全局基准；仅当尺度比失控且无共享边界时才拆分组，组内仍共享基准。

### 步骤 3：法线/切线预算

- 默认 8 位 snorm（分量上界 `1/256`）；要求更高时用 12 位八面体（实测包络 1e−3，`tests.cpp:933`）。
- 含强非均匀缩放的资产，预算要为 inverse-transpose 条件数留余量，并在目标引擎实测角度偏差，不要把分量误差当角度误差用。
- 保持 `reindexTolerance` 与 `nrm_bits` 联动（默认已对齐，`mesh.cpp:776`），不要在外部改 bits 却不更新容差。

### 步骤 4：UV 预算

- 默认 12 位；用 `error/sqrt(uvarea)`（`stream.cpp:253`）作为相对纹素误差指标，超 10% 告警线则 `-vt N` 加位或 `-vtf` 转浮点。
- 高平铺纹理：把平铺倍数计入可接受 UV 绝对误差，必要时为该材质单独提高 `tex_bits`（材质级基准已支持，`stream.cpp:263-279`）。
- 纹理变换保持合并进 `KHR_texture_transform`（`write.cpp:197-207`），不要在量化前烘焙进 UV，以免丧失分组精度。

### 步骤 5：极端输入与合并的防护清单

- 零范围/近零 AABB：依赖现有 `scale==0 → rscale=0` 保护（`stream.cpp:483`），但调用方应剔除真正的点退化网格。
- 超大坐标/远原点：优先 `-vpf` 相对误差路径，避免全局 offset 导致有效位数损失。
- 量化后顶点塌缩：若资产对退化三角形敏感，在量化后增加一次外部焊接/退化剔除（管线本身不做，见 §3.3）。
- morph target：对接缝网格核验基态与 delta 重建后的一致性，必要时对 target 提高位数。
- 重打包：避免对已量化资产二次量化（`parsegltf.cpp:647`），尽量从原始 fp32 资产进入管线。

### 步骤 6：预算冲突的决策层级

1. 库层溢出：接受自动抬 exponent 的精度降级（仅共享指数路径）。
2. 资产层告警：按本方案步骤 2/4 加比特或转浮点，**告警不是失败，而是调用方决策点**。
3. 结构冲突：尊重 `canMergeMeshes` 与并查集的拒绝（不同变换空间/材质不强行共基准）。
4. 视觉回归：简化与量化误差独立叠加，最终以资产在目标引擎的屏幕空间测量为准——这是代码上界无法替代的最后一道验收。

---

## 附：关键代码索引

- 量化原语：`src/quantization.cpp:14-127`、`src/meshoptimizer.h:1182-1202`
- 顶点滤波（Oct/Quat/Exp/Color）：`src/vertexfilter.cpp:78-1499`
- gltfpack 量化基准与告警：`gltf/stream.cpp:131-280, 282-330, 470-659`
- 重索引容差与流程顺序：`gltf/mesh.cpp:770-840, 1324-1353`
- 节点/实例/蒙皮反变换：`gltf/write.cpp:1183-1390`、`gltf/node.cpp:56-187`
- 合并门控：`gltf/mesh.cpp:190-284, 481-529`
- 共享指数 DXR 示例：`demo/main.cpp:935-989`
- 边界与误差测试：`demo/tests.cpp:2842-3043`
- 默认比特与 CLI：`gltf/gltfpack.cpp:1303-1307, 1393-1415`
