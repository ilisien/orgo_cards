from playwright.sync_api import sync_playwright

def generate_pdf():
    url = "http://127.0.0.1:5000/"  # Ensure your Flask app is running
    output_file = "output.pdf"

    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.goto(url, wait_until="networkidle")
        page.pdf(path=output_file, format="A4")
        browser.close()

    print(f"PDF saved as {output_file}")

if __name__ == "__main__":
    generate_pdf()