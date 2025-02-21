
def mirror_backs(pagelists):
    pages = []
    for page in pagelists:
        pages.append([page[2],page[1],page[0],page[5],page[4],page[3],page[8],page[7],page[6]])
    return pages