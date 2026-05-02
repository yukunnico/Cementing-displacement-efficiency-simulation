"""
呼102井 CBL测井PDF文本提取与关键词检索工具
=============================================

本脚本专门针对呼102井CBL（水泥胶结测井）评价相关的PDF文件，
提取全文文本并进行规范化处理，检索与固井质量评价相关的关键词命中行。

与 extract_hu102_sources.py 的区别：
- 本脚本聚焦于CBL评价PDF（100413.PDF等），提取全部页面文本
- 增加了PDF文本的规范化处理（去重复字符、修正标点等）
- 关键词列表针对CBL评价场景定制

输出文件：
1. hu102_cbl_pdf_text_extract.json — 每个PDF的提取文本、字符数和关键词命中
2. hu102_cbl_pdf_key_lines.md — 按文件分组的关键词命中索引（Markdown格式）

用途：从CBL测井报告中提取固井质量合格率、测量井段、油气水层段等关键数据，
为模型验证提供现场对比值。
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

# ========== 路径配置 ==========
ROOT = Path(__file__).resolve().parents[2]  # 项目根目录
PDF_DIR = ROOT / "参考文档" / "呼102" / "1004" / "10041"  # CBL评价PDF所在目录
OUT_DIR = Path(__file__).resolve().parent    # 输出目录
OUT_JSON = OUT_DIR / "hu102_cbl_pdf_text_extract.json"  # 输出：文本提取JSON
OUT_MD = OUT_DIR / "hu102_cbl_pdf_key_lines.md"         # 输出：关键行Markdown

# ========== CBL评价专用关键词 ==========
# 针对CBL固井质量评价报告的关键词，与通用资料扫描的关键词侧重点不同
KEYWORDS = (
    "呼102",         # 井号
    "测量井段",       # CBL测量的深度范围
    "固井质量",       # 固井质量评价
    "目的井段",       # 需要封固的目的层段
    "合格率",         # 固井质量合格率（核心对比值）
    "低于",           # 合格率低于某阈值
    "应封固井段",     # 设计要求封固的井段
    "油气水层段",     # 含油气水层的井段（关键评价对象）
    "悬挂器位置",     # 尾管悬挂器深度
    "完钻井深",       # 完钻井深
    "人工井底",       # 人工井底深度
    "水泥密度",       # 水泥浆密度
    "套管程序",       # 套管结构
    "双层套管",       # 双层套管段（影响CBL评价）
    "不做评价",       # 不做评价的井段说明
)


def collapse_duplicate_chars(text: str) -> str:
    """
    消除PDF提取文本中的重复字符。

    PDF文本提取时常出现字符重复现象（如"固固井井"应为"固井"），
    本函数检测连续两个相同的非ASCII字符，仅保留一个。

    处理逻辑：
    - 如果当前字符和下一个字符相同，且为非ASCII字符（中文字符等），
      则保留一个并跳过下一个
    - ASCII字符（英文、数字）不处理，避免误删合法的重复字母

    Args:
        text: 原始PDF提取文本

    Returns:
        去重后的文本

    示例：
        "固固井井质量量" → "固井质量"
        "hello" → "hello"（ASCII不处理）
    """
    out: list[str] = []
    i = 0
    while i < len(text):
        ch = text[i]
        if i + 1 < len(text) and text[i + 1] == ch and not ch.isascii():
            out.append(ch)
            i += 2  # 跳过重复字符
        else:
            out.append(ch)
            i += 1
    return "".join(out)


def normalize(text: str) -> str:
    """
    PDF文本规范化处理。

    处理步骤：
    1. 消除重复字符（collapse_duplicate_chars）
    2. 修正重复的中文标点符号
    3. 压缩多余空白

    Args:
        text: 原始文本

    Returns:
        规范化后的文本
    """
    text = collapse_duplicate_chars(text)
    # 修正重复的中文标点
    text = text.replace("（（", "（").replace("））", "）")
    text = text.replace("：：", "：").replace("，，", "，").replace("。。", "。")
    # 压缩空白
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def extract_pdf(path: Path) -> dict[str, Any]:
    """
    提取单个PDF文件的文本和关键词命中。

    使用pdfplumber提取全部页面的文本（不限制页数），
    对提取文本进行规范化处理后检索关键词。

    Args:
        path: PDF文件路径

    Returns:
        包含以下字段的字典：
        - path: 文件绝对路径
        - file: 文件名
        - raw_char_count: 原始文本字符数
        - normalized_char_count: 规范化后字符数
        - preview: 规范化文本前4000字符预览
        - hits: 关键词命中列表
    """
    import pdfplumber

    # 提取全部页面文本
    raw_parts: list[str] = []
    with pdfplumber.open(path) as pdf:
        for page in pdf.pages:
            raw_parts.append(page.extract_text() or "")
    raw = "\n".join(raw_parts)

    # 逐行规范化处理
    norm_lines = [normalize(line) for line in raw.splitlines()]
    norm_text = "\n".join(line for line in norm_lines if line)

    # 关键词检索
    hits: list[dict[str, str]] = []
    for idx, line in enumerate(norm_text.splitlines(), start=1):
        matched = [kw for kw in KEYWORDS if kw in line]
        if matched:
            hits.append({"line": str(idx), "keywords": ",".join(matched), "text": line})

    return {
        "path": str(path),
        "file": path.name,
        "raw_char_count": len(raw),
        "normalized_char_count": len(norm_text),
        "preview": norm_text[:4000],
        "hits": hits,
    }


def main() -> None:
    """
    主入口：提取CBL评价PDF文本并生成关键词索引。

    流程：
    1. 扫描PDF目录中的所有.PDF文件
    2. 逐个提取文本、规范化和关键词检索
    3. 保存JSON格式的完整提取结果
    4. 生成Markdown格式的关键词命中索引
    """
    # 提取所有PDF文件
    items: list[dict[str, Any]] = [extract_pdf(path) for path in sorted(PDF_DIR.glob("*.PDF"))]

    # 保存JSON
    OUT_JSON.write_text(json.dumps(items, ensure_ascii=False, indent=2), encoding="utf-8")

    # 生成Markdown关键词索引
    lines = ["# 呼102 CBL PDF 关键行", ""]
    for item in items:
        lines.append(f"## {item['file']}")
        lines.append("")
        for hit in item["hits"]:
            lines.append(f"- L{hit['line']} [{hit['keywords']}]: {hit['text']}")
        lines.append("")
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")

    # 输出统计信息
    print(json.dumps({"pdfs": len(items), "json": str(OUT_JSON), "md": str(OUT_MD)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
