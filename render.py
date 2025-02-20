from flask import Flask, render_template, request

app = Flask(__name__)

# Default value for the template
value = "Hello, World!"

@app.route('/')
def index():
    return render_template('index.html', value=value)

if __name__ == "__main__":
    app.run(debug=True)