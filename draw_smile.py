from sys import argv
from rdkit.Chem import AllChem as AllChem

import rdkit.Chem as Chem
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions #Only needed if modifying defaults
from rdkit.Chem import Draw, MolFromSmiles
import json
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from add import addLine
from arrow import arrowedLine
root="./png/"
path="data.json"
condition_flag = 0
font = ImageFont.truetype('/usr/share/fonts/opentype/noto/NotoSerifCJK-Bold.ttc', 20)  # fc-list :lang=zh
mol_high = 300
mol_width = 500
half_mol_high = int(mol_high/2)
half_mol_width = int(mol_width/2)
one_third_mol_width = int(mol_width/3)
quarter_mol_width = int(mol_width/4)
high = 300
width = 500
half_high=int(high/2)
half_width=int(width/2)
quarter_high=int(high/4)
quarter_width=int(width/4)

add_high = 300
add_width = 50

def draw_reactions(Reactant_list,reaction_index=0,part_num = 1):
    if part_num == 1:
        reaction_width = mol_width
    elif part_num == 2:
        reaction_width = half_mol_width
    elif part_num == 3:
        reaction_width = one_third_mol_width
    elif part_num == 4:
        reaction_width = quarter_mol_width

    try:
        mol_reactant = Chem.MolFromSmiles(Reactant_list[reaction_index])
        img_reaction = Draw.MolToImage(mol_reactant, size=(reaction_width, mol_high),
                                       kekulize=True)  # high=300,whith=500
    except:
        img_reaction = Image.fromarray(np.ones((mol_high, reaction_width, 3), dtype='uint8') * 255)
        draw_reaction = ImageDraw.Draw(img_reaction)
        if len(Reactant_list[reaction_index]) < int(10):
            start = int(reaction_width/2)
            align = "center"
        else:
            start = 0
            align = "left"
        draw_reaction.text((start, half_mol_high), Reactant_list[reaction_index], fill="red", font=font, align=align)
    return img_reaction

with open(path,"r") as data_file:

    for index_line,line in enumerate(data_file):

        path = list(json.loads(line).keys())
        values = list(json.loads(line).values())
        for index_value, value in enumerate(values[0]):
            condition = []
            Reactant_list = []
            Product_list = []
            #统计各role数量
            for role in value[0]:
                if role == "Reactants":
                    for index_rea, Reactant in enumerate(value[0]["Reactants"]):
                        Reactant_list.append(Reactant[0])
                elif role == "Product":
                    Product = value[0]["Product"]
                    Product_list.append(Product[0])
                else:
                    if role == 'Catalyst_Reagents':
                        for index_rea, Catalyst_Reagent in enumerate(value[0]["Catalyst_Reagents"]):
                            condition.append(Catalyst_Reagent[0])
                    elif role == "Yield":
                        condition.append("Yield:"+value[0][role][0][0])
                    else:
                        condition.append(value[0][role][0][0])

            #各部分画图
            #反应物
            add_img= Image.new('RGB', (add_width, add_high), (255, 255, 255))#
            add_img = addLine(add_img, int(add_width/2),int(add_high/2))
            if len(Reactant_list) == 0:
                break
            elif len(Reactant_list) == 1:
                img_reaction = draw_reactions(Reactant_list, reaction_index = 0,part_num = 1)
            elif len(Reactant_list) == 2:
                 part_one = draw_reactions(Reactant_list, reaction_index = 0,part_num = 2)
                 part_two = draw_reactions(Reactant_list,reaction_index = 1,part_num = 2)
                 img_reaction = np.concatenate((part_one, add_img,part_two), axis = 1)
            elif len(Reactant_list) == 3:
                 part_one = draw_reactions(Reactant_list, reaction_index = 0,part_num = 3)
                 part_two = draw_reactions(Reactant_list,reaction_index = 1,part_num = 3)
                 part_three = draw_reactions(Reactant_list, reaction_index = 2,part_num=3)
                 img_reaction = np.concatenate((part_one,add_img, part_two,add_img,part_three), axis = 1)
            elif len(Reactant_list) == 4:
                 part_one = draw_reactions(Reactant_list, reaction_index = 0,part_num = 4)
                 part_two = draw_reactions(Reactant_list,reaction_index = 1,part_num = 4)
                 part_three = draw_reactions(Reactant_list, reaction_index = 2,part_num = 4)
                 part_four = draw_reactions(Reactant_list, reaction_index = 3, part_num = 4)
                 img_reaction = np.concatenate((part_one, add_img, part_two,add_img, part_three,
                                                add_img ,part_four), axis = 1)


            #产物
            if len(Product_list) == 1:
                try:
                    mol_product = Chem.MolFromSmiles(Product_list[0])
                    img_product = Draw.MolToImage(mol_product, size=(mol_width, mol_high), kekulize=True)
                except:
                    img_product = Image.fromarray(np.ones((high, width, 3), dtype='uint8') * 255)
                    draw_product = ImageDraw.Draw(img_product)
                    draw_product.text((half_width, half_high), Product_list[0], fill="red", font=font, align="left")
            # 条件
            if condition !=[]:

                img_condition = Image.fromarray(np.ones((half_high, width, 3), dtype='uint8') * 255)
                draw_condition = ImageDraw.Draw(img_condition)
                con_txt="/".join(condition)
                if len(con_txt) > int(40):
                    con_txt=con_txt[0:20]+"\n"+con_txt[20:-1]
                    start = quarter_width
                    align = "left"
                elif len(con_txt) > int(20):
                    start = 0
                    align = "left"
                draw_condition.text((start, quarter_high), con_txt, fill="red", font=font, align=align)
                condition_flag = 1

            if condition_flag == 1:
                img_arrow = Image.fromarray(np.ones((half_high, width, 3), dtype='uint8') * 255)
                arrow_start,arrow_end = (0, quarter_high), (width, quarter_high)
                img_arrow = arrowedLine(img_arrow, arrow_start, arrow_end, width=1, color=(255, 0, 0))
                img_pad=np.concatenate((img_condition, img_arrow), axis = 0)
                condition_flag = 0
            else:
                img_arrow = Image.fromarray(np.ones((high, width, 3), dtype='uint8') * 255)
                arrow_start,arrow_end = (0, half_high),(width, half_high)
                img_pad = arrowedLine(img_arrow, arrow_start, arrow_end, width=1, color=(255, 0, 0))


            img_rxn= np.concatenate((img_reaction, img_pad, img_product), axis=1)
            image = Image.fromarray(img_rxn)
            image.save(root+ path[0].split("/")[-1].replace(".pdf","") +str(index_line)+ ".png")


