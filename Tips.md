# Seurat scRNA-seq 분석 메모

## Seurat Object 기본 구조

Seurat object는 메모리 절약을 위해 **Logical map** 형태로 데이터를 저장하므로, object를 직접 뜯어서는 내부 정보를 확인할 수 없다.

- `rownames(오브젝트명)` → **유전자 목록** 반환
- `colnames(오브젝트명)` → **Cell Barcode** 반환

`rownames`는 R의 기본 함수이지만, Seurat object에 대해 호출하면 Seurat 라이브러리가 내부적으로 다르게 처리한다. 이를 S4 Method라 한다.

### 주요 메타데이터

- **nFeature_RNA**: 특정 셀에서 발현이 검출된 유전자 종류 수
- **nCount_RNA**: 특정 셀의 총 RNA 분자 수 (UMI로 보정됨)

> VlnPlot에서 dot이 너무 많아 잘 보이지 않을 때는 `pt.size = 0`을 설정하면 된다.

---

## Normalization

### Log Normalization

```r
NormalizeData(데이터명, normalization.method = "LogNormalize", scale.factor = 10000)
```

위가 `NormalizeData`의 기본값이며, Seurat의 표준 정규화 방법이다.

유전자 간 발현량 차이는 매우 극적이기 때문에, 원시 값을 그대로 사용하면 고발현 유전자에 나머지 유전자가 묻힐 수 있다. 이를 보정하기 위해 로그 변환을 적용하고, 값이 지나치게 작아지는 것을 방지하기 위해 10⁴를 곱한 뒤 로그를 씌운다.

### Log Normalization의 한계

유전자 발현이 log-normal distribution을 따른다는 생물학적 근거는 사실 없다. Microarray 시대부터 관습적으로 사용되어 온 방법이다. 실제 scRNA-seq 데이터는 negative binomial distribution에 더 가깝고, 이를 반영한 **SCTransform** 같은 대안이 존재한다. 다만 SCTransform과 Log transform 간의 실질적 차이는 크지 않다고 알려져 있다.

---

## Scaling

### 기본 Scaling

```r
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

`features = all.genes`를 명시하는 이유: 기본값은 VST 분석으로 선정된 variable features 2,000개에 대해서만 scaling을 수행하지만, 이후 Heatmap 등을 그릴 때 다른 유전자도 필요할 수 있으므로 전체 유전자에 대해 scaling하는 것이 일반적이다.

### Regression (공변량 제거)

```r
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
```

variable features에 대해 scaling을 수행하면서, `percent.mt`에 대해 선형 회귀를 적용하여 미토콘드리아 비율에 의한 발현량 차이를 제거한다.

$$x_{ij} = \beta_0 + \beta_1 \cdot \text{percent.mt}_j + \varepsilon_{ij}$$

잔차 $\varepsilon_{ij}$만 남겨 `percent.mt`로 설명되는 변동을 제거하는 원리이며, Cell cycle regression도 동일한 방식이다.

### Regression 적용 시 주의사항

- `percent.mt`에 따른 발현량 변화가 선형적이지 않을 수 있으므로, 최근에는 SCTransform을 권장한다.
- `percent.mt`가 높은 셀은 죽어가는 셀이지만, 이 자체가 생물학적 의미를 가질 수 있다 (예: 항암제 처리 후 데이터).
- Cell cycle 역시 proliferation과 직결되므로, regression을 적용하면 중요한 정보를 잃을 수 있다.
- **Regression은 clustering이 잘 되지 않을 때만 적용하는 것이 좋다.**

---

## PCA (Principal Component Analysis)

scRNA-seq 데이터에는 매우 다양한 종류의 세포가 포함되어 있어, PC1·PC2만으로 2차원 축소를 해서는 의미 있는 클러스터링을 기대하기 어렵다. PCA의 주요 목적은 **PC를 분리하여 UMAP 등 후속 차원 축소 알고리즘의 입력값으로 사용**하는 것이다.

### 왜 PCA를 먼저 하는가

2,000개 variable features를 직접 UMAP에 넣으면 연산량이 과도하다. PCA를 통해 weights의 집합인 PC로 분해한 뒤, 상위 PC(세포 유형 간 주요 변동)와 하위 PC(미세한 차이 또는 noise)를 구분하여 ElbowPlot 등으로 적절한 PC 수를 선택한다.

### 차원 축소 흐름

```
수만 개 유전자 → VST → 2,000개 variable features → PCA → ~20개 미만 PC
```

PCA로 차원을 크게 줄임으로써 **차원의 저주**를 회피하고, 후속 클러스터링에서 유클리드 거리를 효과적으로 활용할 수 있다.

또 PC를 GSEA에 입력하면 각 PC가 어떤 생물학적 차이를 설명하는지 파악할 수 있다.

---

## Clustering

### 주요 함수

```r
FindNeighbors(데이터명, dims = 1:n)   # n개의 PC로 SNN 그래프 구축
FindClusters(데이터명, resolution = r) # 그래프를 r 해상도로 분할
RunUMAP(데이터명, dims = 1:n)          # n개의 PC로 2차원 시각화
```

### RunUMAP의 동작 방식

`RunUMAP`은 단순히 `FindClusters`의 결과를 시각화하는 것이 아니라, 자체적으로 neighbor graph를 구축하고 2차원 거리로 표현한다. 따라서 `RunUMAP`에서도 사용할 PC 수를 별도로 지정해야 한다.

`RunUMAP`과 `FindNeighbors`의 차원 수가 달라도 분석 자체에는 지장이 없지만, UMAP에서 시각적으로 보이는 클러스터와 데이터에 라벨링된 클러스터가 크게 다를 수 있다.

- `CellSelector` 함수로 UMAP 기준 클러스터를 지정할 수도 있으나, UMAP 자체에 왜곡이 있으므로 적절하지 않다.
- 대안으로, PCA 대신 `FindNeighbors`에서 생성한 SNN 그래프를 UMAP 입력으로 사용할 수 있다 (`graph = "RNA_snn"`). 이 경우 UMAP 시각화와 라벨링된 클러스터가 항상 일치한다.

### dims vs resolution
둘다 Cluster를 얼마나 자세히 나눌지 결정하는데 사용할 수 있는 Parameter다.
그러나 의미와 기능에 차이가 있다.

| 파라미터 | 역할 | 영향 |
|---|---|---|
| dims (FindNeighbors) | 셀 간 유사성 계산에 활용할 정보의 양 조절 | 많은 PC → 미묘한 차이까지 포착하지만, 노이즈가 담긴 PC 포함 시 정보 오염 가능. 그래프 자체가 변한다. |
| resolution (FindClusters) | 기존 그래프를 얼마나 잘게 분할할지 결정 | 정보와 그래프는 고정, 커팅 기준만 변경 |

여기서 그래프는 SNN (Shared Nearest Neighbor) graph를 의미한다. 2D에 시각화할 수는 있으나, 셀이 수만 개이고 각 셀이 약 20개의 이웃과 연결되면 edge가 수십만 개에 달하므로 실용적이지 않다.