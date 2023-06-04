import os
from PIL import Image, ImageTk
import tkinter as tk
from tkinter import ttk, filedialog

# set up tkinter window
root = tk.Tk()
root.title("Image Selector")
root.geometry("800x600")

# set up source and destination folders
src_folder = "/Users/hrishijoshi/Desktop/URAP_Research/E4320190426"
dest_folder = "/Users/hrishijoshi/Desktop/URAP_Research/new_images"

# create destination folder if it doesn't exist
if not os.path.exists(dest_folder):
    os.makedirs(dest_folder)

# function to get list of image files in folder
def get_image_files(folder):
    extensions = ["jpg", "jpeg", "png", "gif"]
    files = []
    for file in os.listdir(folder):
        extension = file.split(".")[-1]
        if extension.lower() in extensions:
            files.append(os.path.join(folder, file))
    return files

# create listbox to display image names
listbox = tk.Listbox(root, width=80)
listbox.pack(side=tk.TOP, pady=10)

# function to update listbox with image names
def update_listbox():
    listbox.delete(0, tk.END)
    image_files = get_image_files(src_folder)
    for i, file in enumerate(image_files):
        filename = os.path.basename(file)
        listbox.insert(i, filename)

# function to display image in a new window
def display_image(file):
    image = Image.open(file)
    image.show()

# function to add selected image to destination folder
def add_image(file):
    filename = os.path.basename(file)
    new_file = os.path.join(dest_folder, filename)
    if not os.path.exists(new_file):
        os.makedirs(os.path.dirname(new_file), exist_ok=True)
        with open(new_file, "wb") as f:
            with open(file, "rb") as image:
                f.write(image.read())

# function to handle button click events
def handle_button_click(file):
    if not hasattr(handle_button_click, "display_image_button"):
        handle_button_click.display_image_button = tk.Button(root, text="Display Image", command=lambda: display_image(file))
        handle_button_click.display_image_button.pack(side=tk.LEFT, padx=10)

    if not hasattr(handle_button_click, "add_image_button"):
        handle_button_click.add_image_button = tk.Button(root, text="Add to Destination", command=lambda: add_image(file))
        handle_button_click.add_image_button.pack(side=tk.LEFT, padx=10)

    handle_button_click.display_image_button.configure(command=lambda: display_image(file))
    handle_button_click.add_image_button.configure(command=lambda: add_image(file))

# function to handle label click events
def handle_label_click(event):
    index = listbox.curselection()[0]
    file = get_image_files(src_folder)[index]
    handle_button_click(file)

# bind label click events to handle_label_click function
listbox.bind("<ButtonRelease-1>", handle_label_click)

# create scrollbar for listbox
scrollbar = ttk.Scrollbar(root, orient=tk.VERTICAL, command=listbox.yview)
scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
listbox.config(yscrollcommand=scrollbar.set)

# update listbox with image names
update_listbox()

# allow user to expand the window
root.resizable(True, True)

# run tkinter event loop
root.mainloop()