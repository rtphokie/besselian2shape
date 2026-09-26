from pathlib import Path

import besselian2shape.cache as cache_mod


class _FakeResponse:
    def __init__(self, content: bytes):
        self.raw = _FakeRaw(content)

    def raise_for_status(self) -> None:
        pass


class _FakeRaw:
    def __init__(self, content: bytes):
        self._content = content
        self._served = False

    def read(self, *_args, **_kwargs) -> bytes:
        if self._served:
            return b""
        self._served = True
        return self._content


def test_get_cache_dir_creates_directory(tmp_path, monkeypatch):
    target = tmp_path / "cache_home" / "besselian2shape"
    monkeypatch.setattr(cache_mod, "user_cache_dir", lambda _app: str(target))

    result = cache_mod.get_cache_dir()

    assert result == target
    assert target.is_dir()


def test_download_besselian_csv_downloads_when_missing(tmp_path, monkeypatch):
    monkeypatch.setattr(cache_mod, "get_cache_dir", lambda: tmp_path)

    calls = []

    def fake_get(url, stream, timeout):
        calls.append(url)
        return _FakeResponse(b"year,month,day\n2024,4,8\n")

    monkeypatch.setattr(cache_mod.requests, "get", fake_get)

    path = cache_mod.download_besselian_csv()

    assert path == tmp_path / "eclipse_besselian_from_mysqldump2.csv"
    assert path.read_bytes() == b"year,month,day\n2024,4,8\n"
    assert calls == [cache_mod.BESSELIAN_ELEMENTS_URL]
    assert not path.with_suffix(path.suffix + ".part").exists()


def test_download_besselian_csv_uses_cache_when_present(tmp_path, monkeypatch):
    monkeypatch.setattr(cache_mod, "get_cache_dir", lambda: tmp_path)
    dest = tmp_path / "eclipse_besselian_from_mysqldump2.csv"
    dest.write_bytes(b"cached content")

    def fake_get(*_args, **_kwargs):
        raise AssertionError("should not hit the network when cache exists")

    monkeypatch.setattr(cache_mod.requests, "get", fake_get)

    path = cache_mod.download_besselian_csv()

    assert path == dest
    assert path.read_bytes() == b"cached content"


def test_download_besselian_csv_force_redownloads(tmp_path, monkeypatch):
    monkeypatch.setattr(cache_mod, "get_cache_dir", lambda: tmp_path)
    dest = tmp_path / "eclipse_besselian_from_mysqldump2.csv"
    dest.write_bytes(b"stale content")

    def fake_get(*_args, **_kwargs):
        return _FakeResponse(b"fresh content")

    monkeypatch.setattr(cache_mod.requests, "get", fake_get)

    path = cache_mod.download_besselian_csv(force=True)

    assert path.read_bytes() == b"fresh content"
