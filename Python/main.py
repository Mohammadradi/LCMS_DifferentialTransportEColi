from bs4 import BeautifulSoup as bs
import jinja2
import csv
import os

BATCH_NUMBER = 6

DATA_DIR = f"D:/R_Analysis/output_files/ASKA/Batch_{BATCH_NUMBER}"

INPUT_HTML_FILE = f'Batch{BATCH_NUMBER}_data.html'
OUTPUT_HTML_FILE = f'Batch{BATCH_NUMBER}_Output.html'
PLOT_SUMMARY_FILE = f'allplotdata_{BATCH_NUMBER}.csv'

def get_images(i, soup, data, sign):
    id = f'strain-{i}---{sign}'
    print(id)
    div = soup.find(id=id)
    images = div.find_all('img')
    if len(images) > 0:
        data[i][sign] = list()
        for idx, img in enumerate(images):
            data[i][sign].append({
                'pos': idx+1,
                'src': img['src']
            })


f = open(os.path.join(DATA_DIR, INPUT_HTML_FILE), encoding='utf-8')
soup = bs(f, 'html.parser')

# get the strains
strains = soup.select('#strains code')[1].get_text().split('[1]')[-1].replace('"','').split()
print('Strains', strains)

# get the main div
main_div = soup.find(id='add-the-levels-to-the-plots')
if main_div:
    print('Found main div')

image_dictionary = {}

# get strains div
for j in range(len(strains)):
    i = j + 1
    image_dictionary[i] = dict()
    image_dictionary[i]['neg'] = []
    image_dictionary[i]['pos'] = []
    image_dictionary[i]['name'] = strains[j]
    # get negative images
    get_images(i, soup, image_dictionary, 'neg')
    # get positive images
    get_images(i, soup, image_dictionary, 'pos')


f.close()


# get detailed data
detailed_data = dict()
skip_row = True
with open(os.path.join(DATA_DIR, PLOT_SUMMARY_FILE), encoding='latin-1') as f:
    reader = csv.reader(f)
    for line in reader:
        if skip_row:
            skip_row = False
            continue
        else:
            index, strain, sign, plot, formula, level, name, compound_num = line
            strain, sign, plot = list(map(int, [strain, sign, plot]))
            if sign == 1:
                sign = 'neg'
            else:
                sign = 'pos'
            if strain not in detailed_data:
                detailed_data[strain] = dict()
            if sign not in detailed_data[strain]:
                detailed_data[strain][sign] = dict()
            if plot not in detailed_data[strain][sign]:
                detailed_data[strain][sign][plot] = dict()


            detailed_data[strain][sign][plot] = {
                'formula': formula,
                'level': level,
                'name': name
            }

import json
with open('image.json', 'w') as f:
    json.dump(image_dictionary, f, indent=2)


templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader)
TEMPLATE_FILE = "tabs.html"
template = templateEnv.get_template(TEMPLATE_FILE)
# print(detailed_data)
outputText = template.render(data=image_dictionary, detailed=detailed_data)

with open(os.path.join(DATA_DIR, OUTPUT_HTML_FILE), 'w', encoding='utf-8') as f:
    f.write(outputText)
