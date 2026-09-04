# Off-bin CW leakage in PTA sky maps

当前只研究一个 benchmark：25 pulsars，RA = 6 h、Dec = −45°、fCW = 2.5/T，Earth term only。研究路线见 [WORKFLOW.md](WORKFLOW.md)。

从 [01_MeerKAT_clean_power_hotspot.ipynb](notebooks/01_MeerKAT_clean_power_hotspot.ipynb) 逐 cell 执行。它直接展示 PFOS、Fisher spectrum、正则化决定、clean coefficients 和官方像素转换；已执行至 k = 1–4，并在真实 CW 位置读取 clean / radiometer power。
[00_CW注入_手写原稿.ipynb](notebooks/00_CW注入_手写原稿.ipynb) 保留历史参数和原始输出，不是当前 benchmark；重跑只写缓存。

[03_CW本征频率泄漏.ipynb](notebooks/03_CW本征频率泄漏.ipynb) 从 [GitHub 原稿（466a1ac）](https://github.com/Jiyu-666/Leakage/blob/466a1ac507a0261b9d2889be53876f50b1562da6/demo/03_CW本征频率泄漏.ipynb) 原样恢复，留待后续工作。

## Reproduce

```bash
git submodule update --init --recursive
micromamba create -p .repro-env -f environment.linux-64.lock.yml
micromamba create -p .octave-env -f environment.octave.linux-64.lock.yml
micromamba run -p .repro-env python scripts/prepare_cw_ephemeris.py
micromamba run -p .repro-env python -m ipykernel install --user --name leakage --display-name 'Leakage'
```

在 Jupyter 选择 `Leakage` kernel。Octave 默认从 `.octave-env/bin/octave` 读取；其他安装可用 `LEAKAGE_OCTAVE` 指定。锁文件记录实际执行环境的 Linux x86_64 conda builds（部分依赖要求 x86_64-v3），以及 DEFIANT/pta_replicator 的完整 Git commit；简化环境规格分别是 `environment.yml` 和 `environment.octave.yml`。

PINT 需要完整 DE440 与 BIPM2016 clock data；首次准备需要网络。`prepare_cw_ephemeris.py` 核验完整 DE440，使用 Astropy 的标准缓存。没有换成其他 ephemeris。

运行前必须阅读 [external/COMPATIBILITY.md](external/COMPATIBILITY.md)：指定的 cartography commit 不是开箱即用的现代 Python 包。本项目保存原始 submodule，通过可审计 patch 在 `data/cache/` 生成运行副本。重建和像素转换没有使用旧项目的 sky-map 算法。

```bash
micromamba run -p .repro-env python scripts/validate_cw_notebooks.py test
micromamba run -p .repro-env python scripts/validate_cw_notebooks.py integration
micromamba run -p .repro-env python scripts/validate_cartography.py
```

前两项只验证注入、lineage 与保留的频率响应工具；最后一项验证已生成的 mapping products。当前结果在 `data/products/meerkat_maps/cw_ra6h_dec-45_f2p5T/`，包含共用色标的 [PNG](data/products/meerkat_maps/cw_ra6h_dec-45_f2p5T/clean_power_fbin_01-04.png)、矢量 [PDF](data/products/meerkat_maps/cw_ra6h_dec-45_f2p5T/clean_power_fbin_01-04.pdf) 和真实位置的 [CSV](data/products/meerkat_maps/cw_ra6h_dec-45_f2p5T/power_at_true_cw.csv)。

没有运行 irregular sampling、unequal noise、GWB、pulsar term、分类或显著性分析。Radiometer 只在真实 CW 坐标读取一个标量；不生成 radiometer map，也不运行 S/N 分支。
