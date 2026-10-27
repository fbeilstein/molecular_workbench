"""
SVG Mechanism Normalization and Inkscape Text-to-Path Processor
"""

import os
import re
import shutil
import subprocess
import tempfile


def normalize_svg_text(svg_text: str) -> str:
    """Normalize SVG mechanism text to match Inkscape rendering while preserving arrow strokes and markers."""
    if not svg_text or not svg_text.strip():
        return svg_text

    if "<text" in svg_text:
        inkscape_bin = shutil.which("inkscape")
        if inkscape_bin:
            in_path = None
            out_path = None
            try:
                with tempfile.NamedTemporaryFile("w", suffix=".svg", delete=False) as in_f:
                    in_f.write(svg_text)
                    in_path = in_f.name
                out_path = in_path + ".out.svg"
                cmd = [inkscape_bin, '--actions=select-by-element:text;object-to-path', in_path, '-o', out_path]
                res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=15)
                if res.returncode == 0 and os.path.exists(out_path):
                    with open(out_path, "r") as f:
                        converted = f.read()
                    if "<svg" in converted:
                        return converted
            except Exception as e:
                print(f"[Step7] Inkscape text-to-path error: {e}")
            finally:
                for p in (in_path, out_path):
                    if p and os.path.exists(p):
                        try:
                            os.remove(p)
                        except Exception:
                            pass

        # Fallback regex sanitization if inkscape is not available: only modify <text> elements
        def fix_style(m):
            style = m.group(1)
            style = re.sub(r'stroke:\s*#[0-9a-fA-F]+;?', 'stroke:none;', style)
            if 'font-size' not in style:
                style += ';font-size:10.5px;'
            if 'font-family' not in style:
                style += ';font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;'
            return f'style="{style}"'

        svg_text = re.sub(r'<text\b([^>]*?)style="([^"]*?)"', lambda m: f'<text{m.group(1)}{fix_style(m)}', svg_text)

    return svg_text
