"""Rules: shipped defaults plus an editable user override; evaluate one rule on one metric."""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any

import yaml

DEFAULT_RULES = Path(__file__).with_name("rules.yaml")
USER_RULES = Path.home() / ".config" / "trajinspector" / "rules.yaml"

KINDS = ("max", "min", "range", "not_in", "bool", "manual", "info")
SEVERITIES = ("fail", "warn")
EDITABLE = ("enabled", "value", "severity")   # the only fields the settings table may change

PASS, WARN, FAIL, NA, MANUAL, INFO = "PASS", "WARN", "FAIL", "N/A", "MANUAL", "INFO"


class RulesError(ValueError):
    """The rules file is malformed."""


@dataclass(frozen=True)
class Rule:
    id: str
    gate: int
    check: str
    kind: str
    value: Any
    severity: str
    source: str
    source_url: str | None
    tested_on: str
    rule_of_thumb: bool
    enabled: bool
    note: str = ""
    extra: dict = field(default_factory=dict)   # kind-specific optional keys: warn_in, value_hmr


@dataclass(frozen=True)
class RuleSet:
    version: int
    params: dict
    rules: tuple[Rule, ...]

    def get(self, rule_id: str) -> Rule:
        for rule in self.rules:
            if rule.id == rule_id:
                return rule
        raise KeyError(f"no rule with id '{rule_id}'")


_CORE = {"id", "gate", "check", "kind", "value", "severity", "source", "source_url", "tested_on",
         "rule_of_thumb", "enabled", "note"}


def _parse_rule(raw: dict) -> Rule:
    missing = _CORE - {"note", "source_url"} - raw.keys()
    if missing:
        raise RulesError(f"rule {raw.get('id', '?')}: missing fields {sorted(missing)}")
    if raw["kind"] not in KINDS:
        raise RulesError(f"rule {raw['id']}: unknown kind '{raw['kind']}' (allowed: {KINDS})")
    if raw["severity"] not in SEVERITIES:
        raise RulesError(f"rule {raw['id']}: severity must be one of {SEVERITIES}")
    core = {k: raw.get(k) for k in _CORE if k in raw}
    core.setdefault("source_url", None)
    return Rule(**core, extra={k: v for k, v in raw.items() if k not in _CORE})


def _load_file(path: Path) -> RuleSet:
    try:
        data = yaml.safe_load(path.read_text())
    except (OSError, yaml.YAMLError) as exc:
        raise RulesError(f"cannot read rules file {path}: {exc}") from exc
    rules = tuple(_parse_rule(r) for r in data["rules"])
    ids = [r.id for r in rules]
    if len(ids) != len(set(ids)):
        raise RulesError(f"{path}: duplicate rule ids")
    return RuleSet(int(data["version"]), dict(data["params"]), rules)


def load_rules(user_path: Path | None = USER_RULES) -> RuleSet:
    """Shipped defaults, with enabled/value/severity and params overridden by the user file when it exists.

    Rules missing from the user file keep their shipped values, so a new release can add rules.
    """
    base = _load_file(DEFAULT_RULES)
    if user_path is None or not Path(user_path).exists():
        return base
    user = _load_file(Path(user_path))
    by_id = {r.id: r for r in user.rules}
    merged = tuple(replace(r, **{k: getattr(by_id[r.id], k) for k in EDITABLE}) if r.id in by_id else r
                   for r in base.rules)
    return RuleSet(base.version, {**base.params, **user.params}, merged)


def save_rules(rules: RuleSet, path: Path = USER_RULES) -> None:
    """Write the full rule set so the user file stays readable on its own."""
    def dump(rule: Rule) -> dict:
        d = asdict(rule)
        return {**{k: v for k, v in d.items() if k != "extra"}, **d["extra"]}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump({"version": rules.version, "params": rules.params,
                                    "rules": [dump(r) for r in rules.rules]}, sort_keys=False, allow_unicode=True))


def _violates(rule: Rule, metric: Any) -> bool:
    value = rule.value
    if rule.kind == "max":
        return metric > value
    if rule.kind == "min":
        return metric < value
    if rule.kind == "range":
        return not value[0] <= metric <= value[1]
    if rule.kind == "not_in":
        return metric in value
    return bool(metric) != bool(value)   # bool


def evaluate(rule: Rule, metric: Any) -> str:
    """Status of one metric against one rule. A missing metric must be turned into N/A by the caller."""
    if not rule.enabled:
        return NA
    if rule.kind == "manual":
        return MANUAL
    if rule.kind == "info":
        return INFO
    if metric is None or (isinstance(metric, float) and math.isnan(metric)):
        return NA
    if _violates(rule, metric):
        return rule.severity.upper()
    if rule.kind == "not_in" and metric in rule.extra.get("warn_in", ()):
        return WARN
    return PASS
