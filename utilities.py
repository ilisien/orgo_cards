import requests, os, time
from PIL import Image
from io import BytesIO

def get_placeholder(width=350, height=250):
    """Downloads an image from Picsum and saves it to a file."""
    url = f"https://picsum.photos/{width}/{height}"
    response = requests.get(url)
    
    return Image.open(BytesIO(response.content))

def assemble_card_filesystem(card_title,front_img,back_img):
    directory = f"static/images/{card_title}"
    os.makedirs(directory, exist_ok=True)
    front_path = os.path.join(directory,"front.png")
    back_path = os.path.join(directory,"back.png")
    front_img.save(front_path)
    back_img.save(back_path)

def test_filesystem(card_titles):
    for ct in card_titles:
        assemble_card_filesystem(ct,get_placeholder(),get_placeholder())
        time.sleep(1) #don't overwhelm picsum api

def get_image_paths(card_title,return_which=None):
    # this probably shouldn't be so hardcoded but idk how to fix that atm
    if card_title == "":
        return ""
    elif return_which == None:
        return f"images/{card_title}/front.png", f"images/{card_title}/back.png"
    elif return_which == "back":
        return f"images/{card_title}/back.png"
    elif return_which == "front":
        return f"images/{card_title}/front.png"

def mirror_backs(pagelists):
    pages = []
    for page in pagelists:
        pages.append([page[2],page[1],page[0],page[5],page[4],page[3],page[8],page[7],page[6]])
    return pages

def generate_pagelists(cards):
    while len(cards)%9 != 0:
        cards.append("")
    title_list = [[cards[page*9+card] for card in range(0,9)] for page in range(0,len(cards)//9)]
    fronts = [[get_image_paths(title,"front") for title in page] for page in title_list]
    backs = mirror_backs([[get_image_paths(title,"back") for title in page] for page in title_list])
    return title_list, fronts, backs