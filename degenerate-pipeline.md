# 退化三角形在 meshoptimizer 处理管线中的命运：逐阶段差异来源分析

本文基于仓库实际代码（`src/`、`gltf/`）与可复现的最小实验（链接 `src/*.cpp` 的独立探针，见第 3 节），追踪退化三角形在“合并顶点 → 过滤索引 → 简化 → 空间分簇”各阶段的保留/移除行为，区分三类输入，并解释交换阶段顺序后差异的根因与边界。

---

## 1. 三类“退化”输入与各阶段的判定语义

先把三类容易被混为一谈的输入严格分开。它们在索引层面、位置层面、面积层面的判定条件不同，因此被不同阶段捕获。

| 类别 | 判定条件 | 几何语义 | 典型来源 |
|---|---|---|---|
| **D1 重复索引（拓扑退化）** | 一个三角形内 `a==b || a==c || b==c` | 顶点数 < 3，连索引拓扑都不构成面 | 条带转换残留、手工拼 index buffer |
| **D2 同位置不同属性（接缝楔）** | 三个索引互不相同，但其中两个顶点**位置 bit 级相同**、属性不同 | 面积为零，但顶点键（vertex key）不同，承载 UV/法线接缝语义 | 硬边、UV 缝合线、量化后位置重合的两个楔（wedge） |
| **D3 共线零面积（几何退化）** | 三个索引不同、三个位置也不同，但三点共线 | 叉积为零、有向面积为零，拓扑完全合法 | 数值误差、tessellation 边界、简化过程中相邻顶点被挪到一条线上 |

另需区分两种**非退化**但常被一起讨论的面：

- **重复面（duplicate）**：规范化（保持绕序的循环移位）后顶点键三元组完全相同。
- **反向面（reversed winding）**：同一组顶点、绕序相反。它表达的是另一个有向表面（双面渲染或折叠），**不是冗余**。

### 1.1 阶段一：顶点合并（remap）

涉及两条代码路径，语义不同，是 D2 类差异的第一个分叉点。

**(a) 全顶点键合并** — `meshopt_generateVertexRemap` / `...RemapMulti` / `...RemapCustom`（`src/indexgenerator.cpp:260`）：

- 哈希与相等性比较的对象是**整个顶点字节**（`VertexHasher`，`indexgenerator.cpp:46`）或调用方指定的若干 `meshopt_Stream`。
- gltfpack 的 `reindexMesh`（`gltf/mesh.cpp:792`）走 `generateVertexRemapCustom`：位置必须 bit 级相等，其余属性在 `reindexTolerance`（`mesh.cpp:770`，量化开启时法线 `1/2^nrm_bits`、UV `0.5/2^tex_bits` 等）容差内才视为相等。
- 结论：**D1 原样通过**（remap 只改编号，`(0,1,1)` 仍是 `(r0,r1,r1)`）；**D2 被保留为两个不同顶点**——属性不同导致顶点键不同，接缝楔不会被合并，这正是要保护的语义；D3 也原样通过。

**(b) 纯位置合并** — 简化器内部的 `buildPositionRemap`（`src/simplifier.cpp:201`）与公开的 `meshopt_generatePositionRemap`：

- 只比较前 12 字节位置（`PositionHasher` / `VertexCustomHasher`，相等性见 `indexgenerator.cpp:126`），**完全不看属性**。
- D2 在这一层会被映射到同一个规范位置（探针 Case F：`(0,1,3)` 经位置 remap 变成 `(0,1,0)`），但简化器同时构建 **wedge 环**（`simplifier.cpp:229-239`）把同位置的不同属性顶点串起来，属性接缝信息并没有丢，只是从“不同顶点”变成了“同位置的多个 wedge”。

### 1.2 阶段二：索引过滤

同样有两条路径，判定口径不同。

**(a) gltfpack 的 `filterTriangles`（`gltf/mesh.cpp:842`）**：只做索引相等判断 `a != b && a != c && b != c`。

- 只捕获 **D1**。
- D2（索引不同）、D3（索引不同）、重复面、反向面**全部原样保留**。

**(b) 库函数 `meshopt_filterIndexBuffer[Multi]`（`src/indexgenerator.cpp:339` 的 `filterIndexBuffer`）**，这是语义最完整的过滤器，分两步：

1. 先按调用方给的顶点键（默认全顶点字节；`demo/main.cpp:1497` 用 12 字节位置键）做 remap，再判断 `ra==rb || ...`——因此 **D1 与 D2 都在此被丢弃**（D2 是否丢弃取决于键是否含属性：全顶点键时保留，纯位置键时丢弃，探针 I1 验证了这两种结果）。
2. 把三角形**保持绕序做循环移位规范化**后查哈希表：重复面被丢弃，而**反向面绕序不同、键不同，被保留**（`indexgenerator.cpp:364-383`，头文件注释 `meshoptimizer.h:120` 明确说明 “duplicate triangles with opposite windings are preserved”）。

注意：gltfpack 默认管线（`processMesh`，`gltf/mesh.cpp:1324`）用的是 (a)，**并不调用** `meshopt_filterIndexBuffer`；后者只在 `demo/` 与直接使用库 API 的程序里出现。这意味着经 gltfpack 处理的网格，D2/D3 与重复面都会流入简化器。

### 1.3 阶段三：简化

三条简化路径对退化面的处理方式不同。

**(a) 边折叠简化 `meshopt_simplify[WithAttributes/WithUpdate]`（`src/simplifier.cpp`）**

- 入口先建位置 remap + wedge 表 + 邻接表。退化面会在邻接里产生**自环边**，`classifyVertices` 专门处理：自环会把顶点强制标成 `Kind_Locked`（`simplifier.cpp:390-397`，注释说明否则自环会“闭合”开放边、把顶点误判为流形）。
- 折叠过程中的 `hasTriangleFlips` 跳过退化邻接面（`a == b` 时 continue，`simplifier.cpp:1266`）。
- 真正的清除点在每轮折叠后的 `remapIndexBuffer`（`simplifier.cpp:1832`）：注释写明 “collapse zero area triangles **even if they are not topologically degenerate**”，即基于**位置 remap 后**的 `r0/r1/r2` 比较，D1、D2、D3 只要位置塌缩就会被清掉。
- **关键时序约束（实测确认）**：这个清除函数只在“至少发生过一次边折叠”的轮次后才会执行。若网格因拓扑/误差限制“卡住”（零次折叠），退化面会一直留在输出里。探针 B1/B2/B3：一个好三角形 + 一个 D2 面，即使目标压到 1 个三角形，输出仍保留 2 个面；探针 J（TRACE 版）中存在合法折叠时，`pass 0: removed 1 triangles`，D2 面在该轮 `remapIndexBuffer` 中被清除。
- 属性接缝的保护：`simplifyWithAttributes` 用属性二次误差度量约束折叠；`simplifyUvSplit`（`gltf/mesh.cpp:992`）在简化前把 UV 翻转接缝的顶点**克隆**成独立楔，避免折叠破坏接缝。因此 D2 在“有合法折叠”时通常不是被无差别删掉，而是随其相邻结构一起，在零误差的折叠轮被顺带清除。

**(b) 网格聚类简化 `meshopt_simplifySloppy`（`simplifier.cpp:2646`）**

- 把空间划分为网格单元，同单元顶点取同一 `vertex_id`，再由 `countTriangles`（`simplifier.cpp:2103`）用纯索引比较 `(id0!=id1)&(id0!=id2)&(id1!=id2)` 计数/过滤。
- 它**不看属性、不看绕序语义**：D2 因位置同格而塌缩被删（探针 B-sloppy：2 面 → 1 面，且幸存的是 `(0,1,2)`，接缝楔 v3 所在面消失）；D3 若三个点没落入同一格则会**存活**（探针 C-sloppy：共线面间距 1，网格不够细时 2 面全保留）。头文件与 gltfpack 注释都指出 sloppy “breaks attribute discontinuities”（`gltf/mesh.cpp:1125`）。

**(c) 组件剪除 `meshopt_simplifyPrune`（gltfpack 默认开启，`mesh.cpp:1090`）**

- 按连通分量的包围误差整簇删除，与退化无关，但会让“退化面所在的孤立分量”整体消失——这是分量级而非三角形级的移除，诊断时不要误归为退化过滤。

### 1.4 阶段四：缓存优化、空间排序与 meshlet 分簇

这一阶段的所有操作都是**保面、保绕序的重排/压缩**，不做几何判定：

- `meshopt_optimizeVertexCache[Strip]`：只重排三角形顺序。
- `meshopt_optimizeVertexFetchRemap`（`src/vfetchoptimizer.cpp:7`）：按引用顺序重编号顶点并压缩未引用顶点，三元组的相等关系与绕序不变。
- `meshopt_spatialSortTriangles`（`src/spatialorder.cpp:253`）：按三角形质心做 Morton 序重排，逐面原样搬运（探针 H1：两面绕序 `(0,1,2)`、`(5,4,3)` 排序后绕序不变）。退化面的质心仍是合法点（含重复索引时是两点均值），照样参与排序。
- `meshopt_buildMeshlets`（`src/clusterizer.cpp:1238`）：每个输入三角形恰好发射一次（`emitted_flags` 保证），分簇前后总面数相等（探针 A/B/C 的 `total_tris == input tris`）。退化面只会产生零面积的锥/包围贡献（`demo/tests.cpp:1307` 的 `clusterBoundsDegenerate` 专门覆盖），不会被剔除。
- 点云路径 `meshopt_spatialClusterPoints` / `spatialSortRemap` 只处理顶点，与三角形退化无关。

### 1.5 一张总表：每类输入在哪个阶段被保留/移除

| 输入 | remap（全键/位置） | gltf `filterTriangles` | `filterIndexBuffer`（全键 / 位置键） | 边折叠简化 | sloppy 聚类 | 空间排序 / meshlet |
|---|---|---|---|---|---|---|
| D1 重复索引 | 通过（仅改号） | **移除** | **移除** | 有折叠轮则移除；零折叠则保留 | 同格即移除 | 保留（若活到此时） |
| D2 同位置异属性 | 全键保留 / 位置键合一 | 保留 | 全键**保留** / 位置键**移除** | 有合法折叠时在零误差轮移除；卡住则保留；接缝受属性误差与 UV 克隆保护 | **移除（破坏接缝）** | 保留 |
| D3 共线零面积 | 通过 | 保留 | 保留（键不同） | 有折叠轮则按位置 remap 移除；卡住则保留 | 网格粗时保留、细到同格才移除 | 保留 |
| 重复面 | 通过 | 保留 | **移除**（循环移位规范化后命中） | 不主动去重（实测同绕序重复对在目标压半时仍 2 面共存，仅影响误差度量） | 可能塌缩为同格 | 保留 |
| 反向面 | 通过 | 保留 | **保留**（绕序不同） | 目标不减少时原绕序保留；强制减量时折叠对无合法边、可能被整体移除（见 2.2） | 不区分绕序 | 保留 |

---

## 2. 交换阶段顺序后的差异：预期变换 vs 语义损失

### 2.1 哪些变化属于“预期变换”（无语义差异）

以下变化只改编号或顺序，三角形集合与有向表面完全一致，诊断时应直接忽略：

1. **remap / fetch 重编号**：顶点号变化是任意但确定的双射，`(a,b,c)` 变成 `(r(a),r(b),r(c))`，相等关系与绕序保持。
2. **vertex cache / spatial sort 重排**：三角形在 buffer 中的次序变化，面的内容与绕序不变（探针 H1）。
3. **meshlet 分簇**：面被分组到不同 meshlet，总数不变。
4. **保持绕序的循环移位**：`filterIndexBuffer` 内部把 `(a,b,c)` 规范成最小顶点打头的等价形式，只用于查重，写出时仍用原始索引（`indexgenerator.cpp:387` 写回的是 `a,b,c` 而非规范化键）。

### 2.2 哪些变化构成真实语义损失

1. **属性接缝损失**：把 D2 在“含属性的顶点键”合并之前用**纯位置键**过滤或用 sloppy 聚类，两个接缝楔塌缩成一个，硬边/UV 缝合消失。探针 G 对比：同一网格，位置键 `filterIndexBuffer` 先过滤后只剩 1 面；让 D2 面带着属性进 `simplifyWithAttributes` 则两面俱在。根因是**过滤键的选择**，不是后续简化。
2. **有向表面改变**：
   - 反向面被当作重复面删除（错误实现才会发生；库的规范化只做循环移位、不做翻转，所以 `filterIndexBuffer` 本身安全，头文件 `meshoptimizer.h:120` 有契约）。
   - 简化器的边折叠受 `hasTriangleFlip`（`simplifier.cpp:1231`，法向点积阈值 `0.25*sqrt(...)` ≈ 75°）约束，正常不会翻面；但要注意一个实测边界：当一对反向面构成折叠、且目标面数要求继续减少时，**每条候选边都会触发其中一个面的 flip 检测而全部被判非法**，于是没有任何合法折叠——若此时存在一个零误差的可塌缩结构，折叠轮的 `remapIndexBuffer` 会把这对反向面一起清掉（探针实测：折叠网格目标 1 面时输出 0 面；目标不减少时两面按原绕序完整保留）。这不是翻面重写，而是“无法保持表面时整面移除”；sloppy 路径则根本不做 flip 检查，剧烈聚类可能直接改变表面朝向。
3. **共线面 D3 的存活/消失依赖网格分辨率**：sloppy 的输出随 `target_error` 推导出的 grid size 变化（探针 C 与 J3），这是算法口径差异而非 bug，但在对比两种管线时会表现为面数不同。

### 2.3 最早差异点（root cause）与最终观察结果的关联

把“顺序交换”的两类典型情形拆开：

**情形 A：过滤提前到 remap 之前（或用位置键替代全顶点键）**

- 最早差异点：过滤阶段。此时 D2 的两个楔还没经过全键 remap 确立其“合法的不同顶点”身份，位置键过滤直接把面删掉（探针 I1：同一份数据，全键过滤输出 2 面、位置键过滤输出 1 面）。
- 下游放大：简化器少了一个邻接面，顶点的 `Kind` 分类（border/seam/manifold）随之改变，可选折叠集合不同，最终顶点编号与面数都会不同。**后续 `optimizeVertexFetch` 的编号差异只是症状，根因在过滤键与过滤时机。**
- 正确次序（gltfpack 现状）：先 `reindexMesh`（含属性的自定义键）确立楔，再 `filterTriangles`（只删 D1），把 D2 的取舍权留给“懂属性”的边折叠简化器。

**情形 B：量化提前到 remap 之前**

- gltfpack 实际次序是 `processMesh`（remap/filter/simplify/optimize）全部完成后，才在写出阶段做 KHR_mesh_quantization（`gltfpack.cpp:478` 的 `prepareQuantizationPosition`，量化发生在 `write.cpp` 编码时），即**量化在最后**。
- 若改成先量化：原本 bit 级不同、相距小于半个量化步长的两个位置被量化成同一个值，于是一个原本合法的三角形在量化后变成 D2。探针 E 复现了两条路径：
  - E1（先合并/过滤、后量化）：过滤时面合法、保留；量化后若不再跑第二遍 remap+filter，该面以“量化后零面积”的形态留在最终网格里。
  - E2（先量化、后合并/过滤）：量化使两位置重合，全键 remap 后该面变成 `(0,0,2)` 形态，被拓扑过滤器删除。
- 两条路径最终面数可能相同（E 中都是 1 面），但**语义不同**：E2 是在属性也参与键比较的前提下“知情删除”，E1 残留的是一个渲染时面积为零、却仍占索引与属性插值的面。最早差异点在量化与 remap 的先后，而不是最终 filter。

**情形 C：简化器“卡住”导致的残留**

- 当退化面旁边没有任何合法折叠时（小网格、边界锁、误差目标已满足），`remapIndexBuffer` 永不执行，D1/D2/D3 全部残留到输出，再被空间排序/meshlet 原样带走。实测对比：探针 B 系列（残留 2 面）vs 探针 J 系列（TRACE 显示 `removed 1 triangles`，输出 3 面）。
- 这种残留不能靠后续重排消除；要稳定清除必须在简化**之前**显式过滤（见第 4 节建议）。

---

## 3. 可区分三类原因的最小输入与逐阶段预期结果

下面的输入均可直接喂给库 API；探针程序结构（独立目录、只链 `src/*.o`）可复现全部结果。记顶点为 `位置 | 属性a`。

### 3.1 三个最小输入

```text
输入 A（D1）：
  v0=(0,0,0|1) v1=(1,0,0|1) v2=(0,1,0|1)
  ib = [0,1,2,  0,1,1]                      // 好面 + 重复索引面

输入 B（D2）：
  v0=(0,0,0|1) v1=(1,0,0|1) v2=(0,1,0|1) v3=(0,0,0|9)   // v3 是 v0 的接缝楔
  ib = [0,1,2,  0,1,3]                      // 好面 + 位置退化面

输入 C（D3）：
  v0=(0,0,0|1) v1=(1,0,0|1) v2=(2,0,0|1) v3=(0,1,0|1)   // v0,v1,v2 共线
  ib = [0,1,3,  0,1,2]                      // 好面 + 共线面

输入 D（重复 vs 反向）：
  v0..v2 为普通三角形
  ib = [0,1,2, 0,1,2, 2,1,0]                // 原面 + 精确重复 + 反向面
```

### 3.2 逐阶段预期（实测结果）

| 阶段 | 输入 A | 输入 B | 输入 C | 输入 D |
|---|---|---|---|---|
| 全键 remap 后 | 2 面，`(0,1,1)` 不变 | 2 面，4 顶点全保留 | 2 面（顺序可能被 remap 改） | 3 面 |
| gltf `filterTriangles` 后 | **1 面**（D1 删） | 2 面 | 2 面 | 3 面 |
| `filterIndexBuffer` 全顶点键 | 1 面 | **2 面**（楔保留） | 2 面 | **2 面**：`(0,1,2)` 与 `(2,1,0)` |
| `filterIndexBuffer` 纯位置键 | 1 面 | **1 面**（D2 删） | 2 面 | 2 面 |
| `simplifyWithAttributes`（目标 1 面，无合法折叠） | 1 面 | **2 面（卡住残留）** | 1 面（D3 在唯一折叠轮被清） | — |
| 同左但网格含其它合法折叠（探针 J） | — | D2 在 `pass 0` 零误差清除 | — | — |
| `simplifySloppy` | 1 面 | **1 面（接缝被破坏）** | 粗网格 2 面 / 细网格 1 面 | — |
| `spatialSortTriangles` / `buildMeshlets` | 面数、绕序不变 | 同左 | 同左 | 同左 |

辅助输入 E（量化顺序）与 F（位置 remap 视图）：

```text
输入 E：v0=(0, .|1) v1=(0.02,.|1) v2=(1,.|1) v3=(0,1|1)，量化步长 1/3
  ib = [0,2,3, 0,1,3]
  先量化：v0、v1 重合 → remap 出 (0,0,2) → 过滤后 1 面（知情删除）
  后量化：过滤时 2 面全合法 → 量化后残留一个零面积面（除非再跑一遍 remap+filter）

输入 F：与 B 同数据，仅过 meshopt_generatePositionRemap
  (0,1,3) → (0,1,0)：楔在位置层塌缩，但 wedge 环仍保留属性身份
```

### 3.3 诊断时如何保持原始三角形的对应关系

退化面在管线中会被重编号、重排、甚至拆没，定位根因时建议：

1. **给每个输入三角形分配稳定 ID**（facet id），在每阶段输出时连同三元组一起记录；只比对索引编号必然误判，因为 remap/fetch 重编号是预期行为。
2. **每阶段落盘三份快照**：index buffer、position 流、完整属性流；判定 D2 必须同时看位置与属性，只看索引会把 D2 误当合法面，只看位置会把接缝楔误当重复顶点。
3. **用“面集合 + 绕序”而非“面序列”比较**：先对每个三角形做保持绕序的循环规范化（与 `filterIndexBuffer` 的 `indexgenerator.cpp:364` 同一规则），再比较集合，可消除 spatial sort/cache 优化带来的次序噪声。
4. **开简化器 TRACE 宏**（`src/simplifier.cpp:10` 的 `#define TRACE 1`）重编：`position remap`、`kinds`、每轮 `removed N triangles` 能直接区分“在折叠轮被 `remapIndexBuffer` 清除”与“零折叠残留”。探针 J 的 TRACE 输出即此用法。
5. **分别用全顶点键与位置键各跑一次 `filterIndexBuffer`**：两次输出的差集恰好是 D2 类（位置相同、属性不同）的面，这是分离 D2 与 D1/D3 最快的判别实验。
6. **对量化顺序问题**：在量化前后各做一次位置哈希（`meshopt_generatePositionRemap`），比较 remap 表的等价类变化，新增的等价类就是量化制造的 D2。

---

## 4. 对现有调用约定的兼容处理建议

### 4.1 责任边界：谁该清除退化面

| 层 | 责任 | 不应承担的责任 |
|---|---|---|
| 应用/导出器（gltfpack、demo） | 在进入简化前，用**与自身语义匹配的顶点键**显式过滤 D1，并决定 D2 是否保留 | 不应依赖简化器“顺手”清除（零折叠时不可靠） |
| `meshopt_filterIndexBuffer[Multi]` | 一次性提供 D1+D2（按键）过滤与重复面去重，且契约保证保留反向面 | 不判定 D3（零面积需要几何计算，超出键比较范畴） |
| 边折叠简化器 | 在正常折叠轮通过 `remapIndexBuffer` 清除 D1/D2/D3，并用 wedge/属性误差/UV 克隆保护接缝 | 不保证在“卡住”时清除任何退化面；不负责重复面去重 |
| sloppy 简化器 | 按空间聚类激进降面 | 明确不保护属性接缝与绕序，调用方需自行承担语义损失 |
| 排序/meshlet/fetch 优化 | 只做保面、保绕序的重排与压缩 | 永不删面，诊断时不要把它们当成差异来源 |

### 4.2 具体建议（保持 API 兼容）

1. **gltfpack 侧**：`filterTriangles`（`gltf/mesh.cpp:842`）目前只挡 D1。若要获得确定性输出（不依赖简化器是否卡住），可在 `reindexMesh` 之后、`simplifyMesh` 之前增加一道 `meshopt_filterIndexBufferMulti`，stream 集合与 `reindexMesh` 使用的属性流**完全一致**——这样 D1 与“键级退化”的 D2 被稳定清除，而接缝楔因属性不同仍保留为独立顶点；重复面被去重、反向面按库契约保留。该函数输出长度 ≤ 输入，可原地/同 buffer 调用，不改变现有 `Mesh` 结构与后续签名。
2. **若业务需要保留 D2 接缝面**：维持现状（仅拓扑过滤），但应在文档中明确“简化器零折叠时退化面会残留”，并在写出前可选地再跑一次只针对 D1 的拓扑过滤；不要用纯位置键过滤，否则等于主动放弃接缝。
3. **量化次序**：坚持 gltfpack 现有的“先几何处理、后量化写出”次序。若调用方必须先量化（如自定义管线），应在量化后**补跑一次全键 remap + 过滤**，避免量化制造的零面积面残留（探针 E1 的第二趟）。
4. **D3（共线）无法在键比较层捕获**：需要面积/叉积判定的场景，应在简化前用应用侧几何预过滤，或信任边折叠简化器在有折叠轮时清除；对“零折叠”小网格可在简化后补一道基于位置 remap 的 `meshopt_generatePositionRemap` + 索引相等检查（等价于 `filterIndexBuffer` 的位置键模式），但要明确这一步会同时删除 D2。
5. **重复面与反向面的语义区别必须保留**：去重只能做保持绕序的循环规范化（库的现有行为），切勿在自定义过滤器里把 `(a,b,c)` 与 `(c,b,a)` 归一化到同一键；前者是冗余，后者是有向表面的另一面。探针 D 与 `demo/tests.cpp:2717` 的 `filterTriangles` 测试可作为回归用例。
6. **sloppy 的使用门槛**：仅在允许破坏属性不连续的场合（gltfpack 中对应 `-sa` 且精确简化器卡住的兜底，`mesh.cpp:1126`）使用；含接缝的网格应继续走 `simplifyWithAttributes`。

---

## 附：关键代码位置索引

- 全顶点键 remap / 自定义容差 remap：`src/indexgenerator.cpp:260`、`gltf/mesh.cpp:814`
- 退化 + 重复 + 绕序感知过滤：`src/indexgenerator.cpp:339-394`（契约注释 `src/meshoptimizer.h:116-120`）
- gltfpack 纯拓扑过滤：`gltf/mesh.cpp:842`；管线总入口 `gltf/mesh.cpp:1324`
- 简化器位置 remap 与 wedge 环：`src/simplifier.cpp:201-240`
- 自环退化面的顶点分类：`src/simplifier.cpp:390-397`
- 折叠翻面阈值：`src/simplifier.cpp:1230-1247`
- 零面积面清除（含非拓扑退化）：`src/simplifier.cpp:1832-1864`
- sloppy 聚类计数：`src/simplifier.cpp:2103-2117`、入口 `:2646`
- gltfpack UV 接缝克隆：`gltf/mesh.cpp:992-1062`
- 空间排序保面重排：`src/spatialorder.cpp:253-305`
- meshlet 分簇不删面：`src/clusterizer.cpp:1230-1234`
- 量化阶段位置（写出期）：`gltf/gltfpack.cpp:478`、`gltf/stream.cpp:131`
- 回归测试：`demo/tests.cpp:1973`（simplifyDegenerate）、`:2717`（filterTriangles）、`:1301`（clusterBoundsDegenerate）、`:1856`（sloppy 全退化输入）
