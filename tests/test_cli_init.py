import subprocess
import sys
import json
from pathlib import Path
import pytest

def test_cli_init_creates_structure(tmp_path):
    # مسار إنشاء المشروع التجريبي داخل مجلد مؤقت
    target_dir = tmp_path / "plantvarfilter_project"

    # تنفيذ أمر init لإنشاء ملفات المشروع الأساسية
    result = subprocess.run([
        sys.executable, "-m", "plantvarfilter", "init", str(target_dir)
    ], capture_output=True, text=True)

    # تحقق من نجاح الأمر
    assert result.returncode == 0

    # تحقق من إنشاء المجلدات والملفات المطلوبة
    assert (target_dir / "input").exists()
    assert (target_dir / "output").exists()
    assert (target_dir / "config.json").exists()

    # تحقق من محتويات ملف الإعدادات config.json
    with open(target_dir / "config.json", "r") as f:
        config = json.load(f)

    assert "vcf" in config
    assert "gff" in config
    assert "traits" in config
    assert config["include_intergenic"] is True
    assert isinstance(config["consequence_types"], list)
    assert config["output_format"] in ["csv", "tsv", "json", "feather", "xlsx"]
