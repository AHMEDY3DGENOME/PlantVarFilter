# ui/ui_components.py
import dearpygui.dearpygui as dpg
from plantvarfilter.ui.ui_theme import hgroup, tooltip, COLOR_NOTE

def path_chip_area(tag_group, initial="No file"):
    dpg.add_group(tag=tag_group)
    dpg.add_text(initial, parent=tag_group)

def set_path_chip(tag_group, display_text):
    if dpg.does_item_exist(tag_group):
        dpg.delete_item(tag_group, children_only=True)
        with dpg.group(horizontal=True, horizontal_spacing=6, parent=tag_group):
            dpg.add_text(display_text)
            dpg.add_button(label="×", width=24,
                           callback=lambda: set_path_chip(tag_group, "No file"))

def primary_button(label, callback, tag=None, enabled=True, width=280, height=40):
    return dpg.add_button(label=label, callback=callback, tag=tag,
                          width=width, height=height, enabled=enabled)

def toolbar(title="plantvarfilter", tagline="GWAS & Variant Toolkit"):
    with hgroup(8):
        dpg.add_text(title, tag="brand_title_toolbar")
        tg = dpg.add_text(tagline, tag="brand_tagline_toolbar")
        tooltip(tg, "Variant calling • QC • PLINK • GWAS • GP", color=COLOR_NOTE)
