from flask import Flask, render_template, request

app = Flask(__name__)

# Default value for the template
card_titles = [["number1","number2","number3","number4","number5","number6","number7","number8","number9"],["number10","number11","number12","number13"]]
filenumbers = [[1,2,3,4,5,6,7,8,9],[10,11,12,13]]

@app.route('/')
def index():
    return render_template('flashcards.html', card_titles=card_titles,filenumbers=filenumbers)

if __name__ == "__main__":
    app.run(debug=True)