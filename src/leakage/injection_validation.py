"""Fixed-model raw and centered PINT checks extracted from the old injection notebook."""
import numpy as np
import astropy.units as u
from pint.residuals import Residuals

def cw_weighted_center(values, weights):
    """使用与 PINT 相同的权重去均值；不修改输入数组。"""
    return values - np.average(values, weights=weights)

def cw_error_stats(error):
    """RMS 包括系统偏置，不能用标准差代替。"""
    return {
        "mean_ns": float(np.mean(error)),
        "std_ns": float(np.std(error)),
        "rms_ns": float(np.sqrt(np.mean(error**2))),
        "max_abs_ns": float(np.max(np.abs(error))),
    }

def diagnose_cw_injection(ideal_psr, injected_psr):
    """固定计时模型下的注入检查，不拟合、不修改 TOA/模型/信号。"""
    if ideal_psr.name != injected_psr.name:
        raise ValueError("ideal 与 injected 脉冲星必须同名")
    if ideal_psr.model.as_parfile(include_info=False) != injected_psr.model.as_parfile(include_info=False):
        raise ValueError("注入前后计时模型不一致，不能只比较 CW 增量")
    if "PhaseOffset" in ideal_psr.model.components:
        raise ValueError("本诊断针对隐式 Offset；显式 PHOFF 需要另行处理")
    if not np.array_equal(ideal_psr.toas.table["index"], injected_psr.toas.table["index"]):
        raise ValueError("注入前后 TOA 的数量或行顺序不一致")
    cw_key = f"{injected_psr.name}_cw"
    if set(injected_psr.added_signals_time) != {cw_key}:
        raise ValueError("需要且只能包含一次 CW 注入")

    # 重新计算，而不是直接复用 psr.residuals 的缓存。
    before = Residuals(ideal_psr.toas, ideal_psr.model,
                       subtract_mean=True, use_weighted_mean=True)
    after = Residuals(injected_psr.toas, injected_psr.model,
                      subtract_mean=True, use_weighted_mean=True)
    sigma = before.get_data_error().to_value(u.ns)
    if not np.array_equal(sigma, after.get_data_error().to_value(u.ns)):
        raise ValueError("注入前后 TOA 权重改变")
    if not np.all(np.isfinite(sigma) & (sigma > 0)):
        raise ValueError("TOA uncertainty 必须为有限正数")
    weights = sigma**-2
    baseline = before.time_resids.to_value(u.ns)
    total = after.time_resids.to_value(u.ns)
    cw_raw = injected_psr.added_signals_time[cw_key].to_value(u.ns)
    if cw_raw.shape != baseline.shape or not np.all(np.isfinite(cw_raw)):
        raise ValueError("CW 数组长度不匹配或包含非有限值")
    cw_projected = cw_weighted_center(cw_raw, weights)
    delta_pint = total - baseline
    difference = delta_pint - cw_projected

    # 独立的 raw 检查保留常数偏置，防止去均值掩盖真实误差。
    before_raw = Residuals(ideal_psr.toas, ideal_psr.model,
                           subtract_mean=False, use_weighted_mean=True)
    after_raw = Residuals(injected_psr.toas, injected_psr.model,
                          subtract_mean=False, use_weighted_mean=True)
    delta_raw = after_raw.time_resids.to_value(u.ns) - before_raw.time_resids.to_value(u.ns)
    raw_difference = delta_raw - cw_raw
    legacy_difference = total - cw_raw
    cw_mean = np.average(cw_raw, weights=weights)
    np.testing.assert_allclose(
        legacy_difference, difference + baseline - cw_mean, rtol=0, atol=1e-10
    )
    return {
        "pulsar": injected_psr.name,
        "delta_pint": delta_pint, "cw_projected": cw_projected,
        "difference": difference, "cw_raw": cw_raw,
        "cw_mean_ns": float(cw_mean),
        "baseline": cw_error_stats(baseline),
        "legacy": cw_error_stats(legacy_difference),
        "raw": cw_error_stats(raw_difference),
        "projected": cw_error_stats(difference),
    }
