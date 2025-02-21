from flask import Flask, render_template, request
from utilities import mirror_backs, generate_pagelists

app = Flask(__name__)

# Default value for the template
raw_titles = ["number1","number2","number3","number4","number5","number6","number7","number8","number9","number10","number11","number12","number13"]

card_titles, filenumbers = generate_pagelists(raw_titles)

print(card_titles, filenumbers)

@app.route('/')
def index():
    return render_template('flashcards.html', card_titles=card_titles,filenumbers=filenumbers)

if __name__ == "__main__":
    app.run(debug=True)