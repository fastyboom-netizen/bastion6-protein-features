# Источники

- [Bastion6](https://doi.org/10.1093/bioinformatics/bty155) — определения признаков и классы CTDC/CTDT.
- [PyPredT6](https://github.com/survivebycoding/PyPredT6/blob/master/pypredt6.py) — пример подсчёта состава последовательностей.
- [propy3](https://github.com/MartinThoma/propy3) — матрицы Schneider–Wrede и Grantham для QSO.
- [Данные](https://bastion6.erc.monash.edu/static/download/T6SE_training_data.zip) — учебный набор T6SE.

DPC считает все соседние пары, включая перекрытия, и делится на N−1.
CTDT учитывает оба направления перехода между классами. Промежуточного округления QSO нет,
поэтому значения могут немного отличаться от API propy3.

Матрицы взяты из propy3: Dongsheng Cao, Yizeng Liang, Martin Thoma. Лицензия комплекта — GPL-2.0, полный текст в LICENSE.
