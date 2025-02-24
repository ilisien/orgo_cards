from flask import Flask, render_template, request
from utilities import mirror_backs, generate_pagelists, test_filesystem
from reactions import parse_reaction_file, Reaction

app = Flask(__name__)

#filename = "rxn_lists/openstax/chapter8.txt"
filename = "rxn_lists/test.txt"
reactions = parse_reaction_file(filename)
html_components = []
for reaction in reactions:
    html_components.extend(reaction.all_html_components())

pagified_html_components = generate_pagelists(html_components)
back_pagified_html_components = mirror_backs(pagified_html_components)

print(pagified_html_components)

@app.route('/')
def index():
    return render_template('flashcards.html', phc=pagified_html_components, bphc=back_pagified_html_components)

if __name__ == "__main__":
    app.run(debug=True)