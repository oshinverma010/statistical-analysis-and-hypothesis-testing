"""Download raw data (OWID CO2 + World Bank indicators/metadata) into data/raw.

Rules
-----
* Raw files are cached: an existing file is NEVER overwritten (delete it by hand to force a refresh).
* data/raw/manifest.json records source URL, download date (UTC), size and SHA-256 of each file.
* World Bank calls are retried with back-off; every retry is printed so it can go in the troubleshooting log.
"""
from __future__ import annotations

import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
MANIFEST = RAW / "manifest.json"

OWID_URL = "https://raw.githubusercontent.com/owid/co2-data/master/owid-co2-data.csv"
WB_BASE = "https://api.worldbank.org/v2"
INDICATORS = {
    "NY.GDP.PCAP.KD": "gdp_pc_const2015usd",
    "NY.GDP.PCAP.PP.KD": "gdp_pc_ppp",
    "SP.POP.TOTL": "population_wb",
    "SP.URB.TOTL.IN.ZS": "urban_pct",
    "NV.IND.TOTL.ZS": "industry_pct_gdp",
}
YEAR_RANGE = "1990:2022"
HEADERS = {"User-Agent": "week2-hypothesis-testing-coursework/1.0"}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_manifest() -> dict:
    return json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {}


def _record(manifest: dict, path: Path, url: str, note: str = "") -> None:
    manifest[path.name] = {
        "url": url,
        "downloaded_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S"),
        "bytes": path.stat().st_size,
        "sha256": _sha256(path),
        "note": note,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2))


def _get(url: str, params: dict | None = None, tries: int = 5, timeout: int = 120) -> requests.Response:
    for attempt in range(1, tries + 1):
        try:
            r = requests.get(url, params=params, headers=HEADERS, timeout=timeout)
            r.raise_for_status()
            return r
        except requests.RequestException as exc:  # network error, HTTP error, timeout
            wait = 2 ** attempt
            print(f"  [retry {attempt}/{tries}] {url} -> {type(exc).__name__}: {exc}; sleeping {wait}s")
            time.sleep(wait)
    raise RuntimeError(f"Failed to fetch {url} after {tries} attempts")


def fetch_owid(manifest: dict) -> Path:
    dest = RAW / "owid-co2-data.csv"
    if dest.exists():
        print(f"cached: {dest.name}")
        return dest
    print("downloading OWID CO2 CSV ...")
    r = _get(OWID_URL)
    dest.write_bytes(r.content)
    _record(manifest, dest, OWID_URL)
    return dest


def _wb_paginated(url: str, params: dict) -> list[dict]:
    """Fetch all pages of a World Bank v2 JSON endpoint; raises if the API returns an error message."""
    rows: list[dict] = []
    page = 1
    while True:
        r = _get(url, {**params, "format": "json", "per_page": 20000, "page": page})
        payload = r.json()
        if not isinstance(payload, list) or len(payload) < 2 or payload[1] is None:
            raise RuntimeError(f"Unexpected World Bank response for {url}: {str(payload)[:300]}")
        meta, data = payload
        rows.extend(data)
        if page >= int(meta["pages"]):
            break
        page += 1
    return rows


def fetch_wb_indicator(code: str, manifest: dict) -> Path:
    dest = RAW / f"wb_{code}.csv"
    if dest.exists():
        print(f"cached: {dest.name}")
        return dest
    print(f"downloading World Bank indicator {code} ...")
    url = f"{WB_BASE}/country/all/indicator/{code}"
    rows = _wb_paginated(url, {"date": YEAR_RANGE})
    df = pd.DataFrame(
        {
            "iso3": [r["countryiso3code"] for r in rows],
            "wb_id": [r["country"]["id"] for r in rows],
            "wb_name": [r["country"]["value"] for r in rows],
            "year": [int(r["date"]) for r in rows],
            "value": [r["value"] for r in rows],
        }
    )
    df.to_csv(dest, index=False)
    _record(manifest, dest, url + f"?date={YEAR_RANGE}", note=f"{len(df)} rows incl. aggregates and nulls")
    return dest


def fetch_wb_countries(manifest: dict) -> Path:
    dest = RAW / "wb_country_metadata.csv"
    if dest.exists():
        print(f"cached: {dest.name}")
        return dest
    print("downloading World Bank country metadata ...")
    url = f"{WB_BASE}/country"
    rows = _wb_paginated(url, {})
    df = pd.DataFrame(
        {
            "wb_id": [r["id"] for r in rows],
            "iso3": [r["id"] for r in rows],
            "iso2": [r["iso2Code"] for r in rows],
            "name": [r["name"] for r in rows],
            "region_id": [r["region"]["id"] for r in rows],
            "region": [r["region"]["value"].strip() for r in rows],
            "income_id": [r["incomeLevel"]["id"] for r in rows],
            "income_group": [r["incomeLevel"]["value"] for r in rows],
        }
    )
    df.to_csv(dest, index=False)
    _record(manifest, dest, url, note=f"{len(df)} rows incl. aggregates")
    return dest


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    manifest = _load_manifest()
    fetch_owid(manifest)
    fetch_wb_countries(manifest)
    for code in INDICATORS:
        fetch_wb_indicator(code, manifest)
    print("\nManifest:")
    for name, info in _load_manifest().items():
        print(f"  {name:32s} {info['downloaded_utc']}  {info['bytes']:>10,d} bytes  sha256={info['sha256'][:12]}...")


if __name__ == "__main__":
    main()
