from pathlib import Path
from PIL import Image
from pypdf import PdfWriter

root = Path(__file__).resolve().parents[1]
qa = root / "qa"
qa.mkdir(exist_ok=True)
book = PdfWriter()
for n in range(1, 6):
    pdf, = root.glob(f"Figure_{n}_*_strictTemplate.pdf")
    book.append(pdf)
    with Image.open(pdf.with_suffix(".png")) as img:
        img.thumbnail((990,1100))
        img.convert("RGB").save(qa / f"proof_{n}.png")
book.write(root / "Main_Figures_1_to_5.pdf")
book.close()

