---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.16.1
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Settings

```python
import os

import pandas as pd
```

# Metadata


# Data

```python
data_dir = '/data1/shahs3/users/chois7/projects/GDAN/DLBCL/data/wgs/breakpoints/bedpe_consensus'
```

```python
samples_path = '/data1/shahs3/users/chois7/projects/GDAN/DLBCL/gdan-sv-consensus/resources/samples.txt'
```

```python
samples = [s.strip() for s in open(samples_path).readlines()]
```

```python
agg = pd.DataFrame()
for sample in samples:
    sample_short = '-'.join(sample.split('-')[:4])
    bedpe_path = f'{data_dir}/{sample_short}.SV_consensus.bedpe'
    assert os.path.exists(bedpe_path), bedpe_path
    bedpe = pd.read_table(bedpe_path).rename(columns={'#chrom1':'chrom1'})
    bedpe['sample'] = sample_short
    agg = pd.concat([agg, bedpe])
```

```python
tra = agg[agg['type']=='TRA']
```

```python
tra = tra[['sample', 'gene1', 'gene2']].drop_duplicates()
```

```python
drivers = ['BCL2', 'BCL6', 'MYC']
samples_with_drivers = {g:set() for g in drivers}
cases_with_drivers = {g:set() for g in drivers}
for driver in drivers:
    for _, row in tra.iterrows():
        sample = row['sample']
        case = '-'.join(sample.split('-')[:2])
        gene1 = row['gene1']
        gene2 = row['gene2']
        if gene1 == driver or gene2 == driver:
            samples_with_drivers[driver].add(sample)
            cases_with_drivers[driver].add(case)
```

```python
for driver in drivers:
    print(driver, len(samples_with_drivers[driver]))
```

```python
for driver in drivers:
    print(driver, len(cases_with_drivers[driver]))
```

```python
driver = 'MYC'
tra[(tra['gene1']==driver) | (tra['gene2']==driver)]
```

## Make summary

```python
cases_path = '../resources/cases.txt'
cases = [s.strip() for s in open(cases_path).readlines()]
```

```python
drivers = ['BCL2', 'BCL6', 'MYC']
data = []
for case in cases:
    field = [case]
    for driver in drivers:
        flag = 0
        if case in cases_with_drivers[driver]:
            flag = 1
        field.append(flag)
    data.append(field)
df = pd.DataFrame(data, columns=['case'] + drivers)
```

```python
df.to_csv('../results/v1/summary/rearrangements/driver_rearrangement_status.csv', index=False)
```

```python

```
