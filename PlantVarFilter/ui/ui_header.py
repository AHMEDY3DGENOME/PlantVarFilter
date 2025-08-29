# ui/ui_header.py
import dearpygui.dearpygui as dpg
from PlantVarFilter.ui.ui_theme import hgroup, header_palette

def build_header():
    palette = header_palette()
    with hgroup(10):
        # simple text header (replace with drawlist/logo if you want)
        dpg.add_text("🌱 PlantVarFilter", color=palette["fg"])
        tagline = dpg.add_text("GWAS & Variant Toolkit", color=palette["accent"])
        dpg.add_spacer(width=20)
        dpg.add_text("by Ahmed Yassin", color=palette["muted"])
