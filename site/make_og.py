"""Generate assets/og.png (1200x630), the Open Graph card.

Run by hand when the card text changes: python site/make_og.py
Needs Pillow (not needed by build.py; the PNG is committed). Uses the font that
ships inside Pillow, so the output does not depend on the machine's fonts.
"""
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

OUT = Path(__file__).resolve().parent / "assets" / "og.png"
W, H = 1200, 630
BG, FG, MUTED, ACCENT = (19, 20, 23), (233, 233, 236), (174, 174, 182), (140, 194, 255)

img = Image.new("RGB", (W, H), BG)
d = ImageDraw.Draw(img)
font = lambda size: ImageFont.load_default(size)


def wrap(text, f, width):
    lines, line = [], ""
    for word in text.split():
        trial = f"{line} {word}".strip()
        if d.textlength(trial, font=f) <= width:
            line = trial
        else:
            lines.append(line)
            line = word
    lines.append(line)
    return lines


d.text((64, 56), "THE GRAVITON LADDER", font=font(30), fill=ACCENT)
y = 118
headline = font(56)
for line in ("Entropy cannot tell electric", "tides from magnetic ones.", "Area can."):
    d.text((64, y), line, font=headline, fill=FG)
    y += 72
sub = font(28)
for line in wrap("An exact second-order mismatch for free gravitons in a flat ball", sub, 800):
    d.text((64, y + 14), line, font=sub, fill=MUTED)
    y += 38

d.text((64, 540), "Rerunnable checks. Published failures. Not peer reviewed.", font=font(26), fill=FG)
d.text((64, 578), "ni1ra.github.io/graviton-ladder", font=font(24), fill=MUTED)

# The ladder: entropy / area side on stationary magnetic data, l = 2, 3, 4.
left, right = 900, 1120
top, bottom = 64, 520
for x in (left, right):
    d.line([(x, top), (x, bottom)], fill=ACCENT, width=8)
rungs = [("l = 2", "3/4", 150), ("l = 3", "4/5", 270), ("l = 4", "5/6", 390)]
for tag, ratio, ry in rungs:
    d.line([(left, ry), (right, ry)], fill=ACCENT, width=8)
    f = font(46)
    d.text(((left + right) / 2, ry - 14), ratio, font=f, fill=FG, anchor="ms")
    d.text(((left + right) / 2, ry + 36), tag, font=font(24), fill=MUTED, anchor="ms")
d.text(((left + right) / 2, 566), "entropy / area side", font=font(24), fill=MUTED, anchor="ms")
d.text(((left + right) / 2, 596), "magnetic tides", font=font(24), fill=MUTED, anchor="ms")

OUT.parent.mkdir(exist_ok=True)
img.save(OUT, optimize=True)
print(f"wrote {OUT} ({OUT.stat().st_size} bytes)")
