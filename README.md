# Белковые признаки

Простой модуль для расчёта AAC, DPC, QSO, CTDC и CTDT по последовательности белка.
В папке лежат код, блокнот и готовые CSV для 138 положительных и 1112 отрицательных примеров Bastion6.

[Открыть в Colab](https://colab.research.google.com/github/fastyboom-netizen/bastion6-protein-features/blob/main/Bastion6_features_Colab.ipynb)

В Colab запустите все ячейки и загрузите архив с данными из задания. В конце скачается ZIP с таблицами.

- `protein_features.py` — функции расчёта признаков.
- `Bastion6_features_Colab.ipynb` — загрузка данных, расчёт и выгрузка.
- `positive_*.csv` и `negative_*.csv` — результаты по классам.
- `examples_*.csv` — по одному примеру каждого класса для каждого блока признаков.
- `all_features.csv` — все результаты в одной таблице.
- `feature_dictionary.csv` — названия и описание столбцов.

Всего получается 562 признака: AAC — 20, DPC — 400, QSO — 100, CTDC — 21, CTDT — 21.
Частоты записаны долями. Для QSO используются `maxlag=30` и `weight=0.1`.

```python
from protein_features import calculate_features

sequence = "ACDEFGHIKLMNPQRSTVWY" * 2
features = calculate_features(sequence)
```

Источники и лицензия указаны в `SOURCES.md` и `LICENSE`.
