# ui/ui_header.py
import dearpygui.dearpygui as dpg
from plantvarfilter.ui.ui_theme import hgroup, header_palette
import os

def _build_header(self, parent, big: bool = False):
    with dpg.group(parent=parent, horizontal=True, horizontal_spacing=10):
        logo_path = os.path.join(os.path.dirname(__file__), "assets", "logo.png")
        if os.path.exists(logo_path):
            w, h, c, data = dpg.load_image(logo_path)
            if not dpg.does_item_exist("plant_logo_tex"):
                with dpg.texture_registry(show=False):
                    dpg.add_static_texture(width=w, height=h, default_value=data, tag="plant_logo_tex")
            dpg.add_image("plant_logo_tex", width=52 if big else 40, height=52 if big else 40)
        else:
            dl = dpg.add_drawlist(width=52 if big else 40, height=52 if big else 40)
            dpg.draw_circle(center=(26 if big else 20, 26 if big else 20),
                            radius=24 if big else 18,
                            color=(76, 175, 110, 255), thickness=2, parent=dl)

        # العنوان
        title = dpg.add_text("plantvarfilter",
             color=(210, 230, 210) if self.night_mode else (30, 45, 35),
             default_font=self._font_title)

        if self._font_title:
            dpg.bind_item_font(title, self._font_title)

        dpg.add_spacer(width=10)
        tagline = dpg.add_text("GWAS & Variant Toolkit",
                               color=(160, 200, 160) if self.night_mode else (40, 90, 40))
        if self._font_subtitle:
            dpg.bind_item_font(tagline, self._font_subtitle)

        dpg.add_spacer(width=14)
        author = dpg.add_text("by Ahmed Yassin",
                              color=(140, 150, 140) if self.night_mode else (50, 70, 55))
        if self._font_subtitle:
            dpg.bind_item_font(author, self._font_subtitle)

