def mirror_backs(pagelists):
    pages = []
    for page in pagelists:
        pages.append([page[2],page[1],page[0],page[5],page[4],page[3],page[8],page[7],page[6]])
    return pages

def generate_pagelists(cards):
    while len(cards)%9 != 0:
        cards.append(-1)
    title_list = [[cards[page*9+card] for card in range(0,9) if cards[page*9+card] != -1] for page in range(0,len(cards)//9)]
    filenumbers = mirror_backs([[str(page*9+card+1) if cards[page*9+card] != -1 else str(-1) for card in range(0,9) ] for page in range(0,len(cards)//9)])
    return title_list, filenumbers
