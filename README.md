# 分子构象几何核查工具 (Molecular Geometry Audit)

使用 Python 标准库（无任何第三方依赖、不依赖化学数据库）读取分子拓扑与
多组构象坐标，计算局部几何测量、整体大小指标和相对参考构象的距离差指标，
定位可疑结构（零长键、共线链等导致的未定义测量）。

## 调用方式

```bash
python -m molecular_geometry \
    --coordinates examples/coordinates.csv \
    --topology  examples/topology.json \
    --config     examples/config.json \
    --output-dir examples/out
```

成功时在输出目录写入 `geometry.csv` 与 `analysis.json`（输出目录不存在会自动
创建）。任何输入校验失败均以非零状态码退出（CLI 退出码 2），且**不会创建或
覆盖**任何报告文件（先完成全部校验与计算，再原子写入）。

## 输入格式

### 拓扑 JSON

```json
{
  "atoms": [
    {"id": "A1", "mass": 1.0},
    {"id": "A2", "mass": 2.0}
  ],
  "bonds": [
    {"atoms": ["A1", "A2"]}
  ]
}
```

- 原子编号为非空字符串；质量必须为正的有限数值。
- 键为无向键，按原子编号对归一化排序；`["A1","A2"]` 与 `["A2","A1"]`
  视为同一条键。工具不推测键，也不推断键级。
- 拒绝：空拓扑、重复原子、非正/非有限质量、自连键、重复键（含端点颠倒）、
  引用未知原子的键。

### 坐标 CSV

首行为表头，必须包含 `conformer`、`atom`、`x`、`y`、`z` 列（表头大小写
不敏感；构象列也接受 `conf` 等常见别名）。

- 每个构象必须恰好覆盖拓扑中的全部原子。
- 拒绝：空表、同一构象内原子重复、原子缺失/未知、坐标非数值或非有限。
- 坐标单位由配置全局指定：`"angstrom"`（默认）或 `"nm"`；纳米值乘以 10
  统一换算为埃。

### 配置 JSON

```json
{"reference_conformer": "ref", "coordinates_unit": "angstrom"}
```

`reference_conformer` 必须是坐标表中已存在的构象。

## 输出

### geometry.csv

按构象编号、再按原子编号序列字典序排序；列：

| 列 | 含义 |
| --- | --- |
| `conformer` | 构象编号 |
| `measurement` | `bond_length` / `bond_angle` / `torsion` |
| `atoms` | 原子编号序列（JSON 数组） |
| `value` | 数值；未定义时为空（**绝不输出 0 或 NaN**） |
| `unit` | `angstrom` 或 `degree` |
| `undefined_reason` | 未定义原因；有值时为空字符串 |

### analysis.json

- `units`：长度=埃，角度=度，并记录输入坐标单位。
- `geometry_conventions`：各项测量与指标的口径说明。
- 每个构象包含：
  - `center_of_mass_angstrom`：**质量**中心（不是几何中心）；
  - `radius_of_gyration_angstrom`：质量加权回转半径，即各原子到质量中心
    距离平方按质量加权平均后的平方根；
  - `measurement_counts`：各类型及合计的有效/未定义测量数量；
  - `bond_length_angstrom`、`bond_angle_degrees`：有效值的最小、最大、均值，
    无有效值时三者均为 `null`；
  - `distance_difference_rms_angstrom` 与 `max_distance_difference`。

## 几何口径

- **键长**：成键两原子的欧氏距离；键长本身允许为零。
- **键角**：共享中心原子的相邻两键，从中心原子分别指向两个端点的外向向量
  的夹角，范围 0–180°；按端点编号排序，只报告一次。
- **扭转角**：由四个不同原子构成的三键链 a-b-c-d，取前两键向量叉积与后
  两键向量叉积（即两个平面的法向）之间的**无符号**夹角，范围 0–180°；
  反向链 d-c-b-a 以原子编号序列较小者为准，只报告一次。
- **未定义项**：零长键使相关键角/扭转角未定义；连续三原子共线使该链扭转
  角未定义。未定义返回 `null` 并写明原因，其他测量照常完成。
- **参考距离差 RMS**：枚举全部无序原子对（编号并列选取），对构象内原子对
  距离与参考构象对应距离之差求均方根；并报告绝对差最大的原子对（并列时按
  编号取较小者）。少于两个原子时返回 `null` 并说明原因。该指标基于内部
  距离，**平移或旋转构象不改变其值**，不直接比较未经对齐的坐标。它**不是**
  手性鉴别手段，也**不是**能量评价。

## 示例（`examples/`）

- `ref`：三个互相垂直的单位键，键角与扭转角均为 90°。
- `linear`：直线结构，键角 180°，扭转角因共线而未定义（null + 原因）。
- `coincident`：两个原子坐标重合，出现零长键；对应键长仍为 0，相关键角与
  扭转角未定义。
- `rotated`：`ref` 的刚体平移；所有内部几何不变，参考距离差 RMS 为 0。

运行：

```bash
python -m molecular_geometry --coordinates examples/coordinates.csv \
    --topology examples/topology.json --config examples/config.json \
    --output-dir examples/out
```

## 测试

```bash
python -m unittest discover -s tests -v
```
