from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
from PIL import Image

# Constants
CARD_WIDTH = 2.5 * 72  # Convert inches to points (1 inch = 72 points)
CARD_HEIGHT = 3.5 * 72
PAGE_WIDTH, PAGE_HEIGHT = letter  # 8.5 x 11 inches
MARGIN = 20  # Margin from page edge
CARDS_PER_ROW = 3
CARDS_PER_COLUMN = 3

flashcards = [
    ("Front 1", "back1.png"),
    ("Front 2", "back2.png"),
    ("Front 3", "back3.png"),
    ("Front 4", "back4.png"),
    ("Front 5", "back5.png"),
    ("Front 6", "back6.png"),
    ("Front 7", "back7.png"),
    ("Front 8", "back8.png"),
    ("Front 9", "back9.png"),
]

def draw_flashcard(c, x, y, front_text=None, front_image=None):
    """ Draws a single flashcard with text or an image """
    c.rect(x, y, CARD_WIDTH, CARD_HEIGHT)  # Card border
    
    if front_text:
        c.setFont("Helvetica", 12)
        c.drawCentredString(x + CARD_WIDTH / 2, y + CARD_HEIGHT / 2, front_text)
    elif front_image:
        img = ImageReader(front_image)
        c.drawImage(img, x, y, CARD_WIDTH, CARD_HEIGHT, preserveAspectRatio=True, anchor='c')

def create_flashcard_pdf(output_filename, flashcards):
    """ Generates a PDF with flashcards arranged on two pages (fronts & backs) """
    c = canvas.Canvas(output_filename, pagesize=letter)
    
    index = 0
    for row in range(CARDS_PER_COLUMN):
        for col in range(CARDS_PER_ROW):
            if index >= len(flashcards):
                break
            x = MARGIN + col * (CARD_WIDTH + 10)
            y = PAGE_HEIGHT - (MARGIN + (row + 1) * (CARD_HEIGHT + 10))
            
            front_text, _ = flashcards[index]
            draw_flashcard(c, x, y, front_text=front_text)
            index += 1

    c.showPage()  # Start a new page for the backs

    index = 0
    for row in range(CARDS_PER_COLUMN):
        for col in range(CARDS_PER_ROW):
            if index >= len(flashcards):
                break
            # Mirror horizontally by reversing the column order
            mirrored_col = (CARDS_PER_ROW - 1 - col)
            x = MARGIN + mirrored_col * (CARD_WIDTH + 10)
            y = PAGE_HEIGHT - (MARGIN + (row + 1) * (CARD_HEIGHT + 10))
            
            _, back_image = flashcards[index]
            draw_flashcard(c, x, y, front_image=back_image)
            index += 1

    c.showPage()
    c.save()

# Generate the PDF
create_flashcard_pdf("flashcards.pdf", flashcards)
