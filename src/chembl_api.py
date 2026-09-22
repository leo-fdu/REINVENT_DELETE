#!/usr/bin/env python3
"""Minimal ChEMBL REST helper with retries, rate limiting, and pagination."""

from __future__ import annotations

import time
from typing import Any, Dict, Iterator, Optional

import requests

CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data"
DEFAULT_TIMEOUT = 90.0
DEFAULT_MAX_RETRIES = 16
DEFAULT_BACKOFF = 3.0
MAX_BACKOFF = 90.0
MIN_REQUEST_INTERVAL = 3.0

_last_request_time = 0.0


def _throttle() -> None:
    """Space out requests to stay clear of ChEMBL rate limits."""
    global _last_request_time
    elapsed = time.time() - _last_request_time
    if elapsed < MIN_REQUEST_INTERVAL:
        time.sleep(MIN_REQUEST_INTERVAL - elapsed)
    _last_request_time = time.time()


def get_json(
    url: str,
    params: Optional[Dict[str, Any]] = None,
    timeout: float = DEFAULT_TIMEOUT,
    max_retries: int = DEFAULT_MAX_RETRIES,
    backoff: float = DEFAULT_BACKOFF,
    session: Optional[requests.Session] = None,
) -> Dict[str, Any]:
    """GET a ChEMBL JSON resource, retrying on 5xx/timeouts with capped backoff."""
    http = session or requests.Session()
    last_error: Optional[Exception] = None
    for attempt in range(max_retries):
        _throttle()
        try:
            response = http.get(url, params=params, timeout=timeout)
        except (requests.Timeout, requests.ConnectionError) as exc:
            last_error = exc
            sleep_seconds = min(backoff * (2**attempt), MAX_BACKOFF)
            print(f"  [retry {attempt + 1}/{max_retries}] {exc} -> sleep {sleep_seconds:.0f}s")
            time.sleep(sleep_seconds)
            continue
        if response.status_code in (429, 500, 502, 503, 504):
            last_error = RuntimeError(f"HTTP {response.status_code}")
            sleep_seconds = min(backoff * (2**attempt), MAX_BACKOFF)
            print(
                f"  [retry {attempt + 1}/{max_retries}] HTTP {response.status_code} "
                f"-> sleep {sleep_seconds:.0f}s"
            )
            time.sleep(sleep_seconds)
            continue
        response.raise_for_status()
        return response.json()
    raise RuntimeError(f"ChEMBL request failed after {max_retries} attempts: {url}") from last_error


def iter_activity_pages(
    target_chembl_id: str,
    page_size: int = 1000,
    start_offset: int = 0,
    session: Optional[requests.Session] = None,
) -> Iterator[Dict[str, Any]]:
    """Yield every activity record for one target, page by page."""
    offset = start_offset
    total_count: Optional[int] = None
    while total_count is None or offset < total_count:
        payload = get_json(
            f"{CHEMBL_API_BASE}/activity.json",
            params={
                "target_chembl_id": target_chembl_id,
                "limit": page_size,
                "offset": offset,
            },
            session=session,
        )
        records = payload.get("activities", [])
        total_count = payload["page_meta"]["total_count"]
        if not records:
            break
        for record in records:
            yield record
        offset += len(records)
        if offset < total_count:
            print(f"    {target_chembl_id}: {offset}/{total_count}")


def fetch_molecule(
    molecule_chembl_id: str,
    session: Optional[requests.Session] = None,
) -> Optional[Dict[str, Any]]:
    """Fetch a single molecule record, returning None on persistent 404/failure."""
    try:
        return get_json(f"{CHEMBL_API_BASE}/molecule/{molecule_chembl_id}.json", session=session)
    except Exception as exc:
        print(f"  [warn] molecule {molecule_chembl_id}: {exc}")
        return None


def fetch_molecules(
    molecule_chembl_ids: list[str],
    batch_size: int = 10,
    session: Optional[requests.Session] = None,
) -> Iterator[Dict[str, Any]]:
    """Yield molecule records in small batches, falling back to single lookups."""
    total = len(molecule_chembl_ids)
    for start in range(0, total, batch_size):
        batch = molecule_chembl_ids[start : start + batch_size]
        payload = None
        if len(batch) > 1:
            try:
                payload = get_json(
                    f"{CHEMBL_API_BASE}/molecule/set/{';'.join(batch)}.json",
                    session=session,
                    max_retries=3,
                )
            except Exception as exc:
                print(
                    f"  [warn] molecule/set batch of {len(batch)} failed ({exc}); "
                    f"falling back to single lookups"
                )
        if payload is not None:
            for record in payload.get("molecules", []):
                yield record
        else:
            for molecule_id in batch:
                record = fetch_molecule(molecule_id, session=session)
                if record is not None:
                    yield record
        print(f"  molecules: {min(start + len(batch), total)}/{total}")
