from flask import Flask, render_template, request
from utilities import mirror_backs, generate_pagelists, test_filesystem

app = Flask(__name__)

# Default value for the template
raw_titles = ["number1","number2","number3","number4","number5","number6","number7","number8","number9","number10","number11","number12","number13"]

card_titles, fronts, backs = generate_pagelists(raw_titles)

print(card_titles, fronts, backs)

@app.route('/')
def index():
    return render_template('flashcards.html', card_titles=card_titles,fronts=fronts,backs=backs)

if __name__ == "__main__":
    app.run(debug=True)