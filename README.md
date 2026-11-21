# Белковые признаки Bastion6

Python-модуль и блокнот Google Colab для расчёта частотных и физико-химических признаков белковых последовательностей. Результаты получены для учебного набора T6SE: **138 положительных и 1112 отрицательных примеров**.

[Открыть блокнот из репозитория в Colab](https://colab.research.google.com/github/fastyboom-netizen/bastion6-protein-features/blob/main/Bastion6_features_Colab.ipynb) · [Выполненная версия на Google Диске](https://colab.research.google.com/drive/15Rzfk9St41wJR2ZC6OLKLz3hKonqCxbI)

Ссылка на Google Диск использует настройки доступа владельца. Для повторного расчёта доступ к ней не требуется: запускаемая версия блокнота находится в этом репозитории без сохранённых выводов. CSV и `manifest.json` взяты из комплекта, созданного при выполнении в Colab.

## Рассчитываемые признаки

| Блок | Размер | Определение |
|---|---:|---|
| AAC | 20 | Частота каждой аминокислоты: `n(a)/N` |
| DPC | 400 | Частота каждой упорядоченной соседней пары: `n(ab)/(N−1)`, с перекрытиями |
| QSO | 100 | По 50 компонент для матриц Schneider–Wrede и Grantham: 20 частотных и 30 порядковых |
| CTDC | 21 | Доли трёх классов каждого из семи физико-химических свойств |
| CTDT | 21 | Частоты межклассовых переходов `1↔2`, `1↔3`, `2↔3` для семи свойств |
| **Всего** | **562** | Признаки каждой последовательности |

CTDC и CTDT используют гидрофобность, нормализованный объём ван-дер-Ваальса, полярность, поляризуемость, заряд, вторичную структуру и доступность растворителю. Классы взяты из таблицы 1 Bastion6 и явно заданы в `GROUPS` модуля.

Для QSO используются `maxlag=30`, `weight=0.1`. Для каждого лага `d` рассчитывается сумма квадратов расстояний между остатками на позициях `i` и `i+d`. Нормирующий знаменатель равен `1 + weight × sum(tau)`. Матрицы включены в модуль, загрузка из сети при расчёте не нужна. Все частоты выражены долями; промежуточного округления нет.

## Запуск в Colab

1. Откройте блокнот по ссылке выше.
2. Выберите **Среда выполнения → Выполнить все**. Достаточно CPU.
3. Загрузите ZIP с [учебными данными Bastion6](https://bastion6.erc.monash.edu/static/download/T6SE_training_data.zip).
4. После выполнения скачайте `bastion6_features_bundle.zip`.

Блокнот создаёт модуль, таблицы по блокам и классам, общую таблицу, словарь признаков и манифест. В новой сессии файлы создаются в `bastion6_run/results/`. В репозитории готовые CSV размещены в корне для прямого скачивания.

Для блокнота нужны NumPy и pandas, уже доступные в Colab. Модуль `protein_features.py` использует только стандартную библиотеку Python. BioPython и propy3 устанавливать не требуется.

## Использование модуля

```python
from protein_features import calculate_features, aac, dpc, qso, ctdc, ctdt

sequence = "ACDEFGHIKLMNPQRSTVWY" * 2
features = calculate_features(sequence, maxlag=30, weight=0.1)
print(len(features))  # 562
print(aac(sequence)["AAC_A"])  # 0.05
```

Для QSO требуется `N > maxlag`. Пустые строки, неизвестные остатки, пропуски и стоп-символы вызывают ошибку. Регистр и пробельные символы нормализуются. Последовательности не обрезаются, остатки не удаляются. Все соседние пары учитываются, включая перекрытия. Встроенный FASTA-парсер поддерживает многострочные последовательности и сохраняет последнюю запись.

## Примеры и полные результаты

В каждом файле `examples_*` приведены один положительный и один отрицательный белок со всеми признаками соответствующего блока.

| Блок | Два примера | Положительные: 138 | Отрицательные: 1112 |
|---|---|---|---|
| AAC | [examples_AAC.csv](examples_AAC.csv) | [positive_AAC.csv](positive_AAC.csv) | [negative_AAC.csv](negative_AAC.csv) |
| DPC | [examples_DPC.csv](examples_DPC.csv) | [positive_DPC.csv](positive_DPC.csv) | [negative_DPC.csv](negative_DPC.csv) |
| QSO | [examples_QSO.csv](examples_QSO.csv) | [positive_QSO.csv](positive_QSO.csv) | [negative_QSO.csv](negative_QSO.csv) |
| CTDC | [examples_CTDC.csv](examples_CTDC.csv) | [positive_CTDC.csv](positive_CTDC.csv) | [negative_CTDC.csv](negative_CTDC.csv) |
| CTDT | [examples_CTDT.csv](examples_CTDT.csv) | [positive_CTDT.csv](positive_CTDT.csv) | [negative_CTDT.csv](negative_CTDT.csv) |

- [all_features.csv](all_features.csv) — все 1250 белков и 562 признака.
- [positive_features.csv](positive_features.csv), [negative_features.csv](negative_features.csv) — все признаки по классам.
- [feature_dictionary.csv](feature_dictionary.csv) — определение каждого столбца.
- [sequences.csv](sequences.csv) — последовательности и исходные заголовки FASTA.
- [manifest.json](manifest.json) — параметры, версии, размеры выборок и SHA-256.

CSV записаны с разделителем `,`, десятичной точкой и кодировкой UTF-8 с BOM. Первые девять столбцов — метаданные; остальные — числовые признаки. Метки проекта: `positive=1`, `negative=0`. `record_key` идентифицирует строку, `sequence_sha256` позволяет сопоставлять последовательности. Автоматическое удаление повторов не выполняется; в использованном наборе повторов последовательностей не обнаружено.

## Проверка результатов

В Colab успешно обработаны все 1250 записей. Выполнены проверки на контрольных последовательностях, включая перекрывающиеся пары `AAAA`, переходы классов заряда и QSO для однородной и неоднородной последовательностей. Проверены размер вектора, конечность и неотрицательность значений, суммы AAC/DPC/CTDC и обоих векторов QSO, границы CTDT и обратное чтение CSV. Параметры окружения сохранены в манифесте: Python 3.13.18, NumPy 2.4.6, pandas 2.2.3.

## Источники и соглашения

1. [Bastion6: a bioinformatics approach for accurate prediction of type VI secreted effectors](https://doi.org/10.1093/bioinformatics/bty155), разделы 2.2.1, 2.2.3 и таблица 1.
2. [PyPredT6: pypredt6.py](https://github.com/survivebycoding/PyPredT6/blob/master/pypredt6.py).
3. [propy3](https://github.com/MartinThoma/propy3) — источник матриц QSO.
4. [Документация QSO](https://propy3.readthedocs.io/en/latest/QuasiSequenceOrder.html).

Расчёт следует определениям Bastion6. В PyPredT6 DPC использует `str.count` без перекрытий и другой знаменатель; эти особенности здесь не воспроизводятся. propy3 округляет промежуточные величины QSO, поэтому побитового совпадения с его API не предполагается. Подробности — в [SOURCES.md](SOURCES.md).

Обучение классификатора и эволюционные признаки не входят в эту работу. Комплект кода и матрицы распространяются под [GPL-2.0](LICENSE). Условия использования исходного набора данных определяются его источником.
