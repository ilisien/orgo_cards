from weasyprint import HTML

HTML("flashcards.html").write_pdf("flashcards.pdf")
