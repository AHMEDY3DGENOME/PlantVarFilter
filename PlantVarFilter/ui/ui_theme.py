# ui/ui_theme.py
import os
import dearpygui.dearpygui as dpg
from contextlib import contextmanager

# ---------------- Design tokens ----------------
FONT_SCALE = 1.10
SPACING = 14
RADIUS = 10
CARD_W, CARD_H = 560, 280

# ---------------- Palette ----------------
COLOR_TEXT_DARK = (234, 238, 241)
COLOR_TEXT_LIGHT = (28, 32, 33)

COLOR_BG_DARK_WINDOW = (22, 25, 26)
COLOR_BG_DARK_CHILD  = (30, 34, 35)
COLOR_BG_DARK_FRAME  = (42, 46, 47)

COLOR_BG_LIGHT_WINDOW = (246, 247, 248)
COLOR_BG_LIGHT_CHILD  = (236, 239, 242)
COLOR_BG_LIGHT_FRAME  = (225, 228, 232)

COLOR_TITLE  = (205, 228, 210)
COLOR_NOTE   = (219, 191, 84)
COLOR_MUTED  = (152, 163, 169)

# Default primary (will be overridden at runtime)
COLOR_PRIMARY   = (46, 125, 50)
COLOR_PRIMARY_H = (67, 160, 71)
COLOR_PRIMARY_A = (27, 94, 32)

# ---------------- Font helpers ----------------
_APP_FONT_TAG = "ui_font_base"
_FONT_REGISTRY_TAG = "font_registry"

def _project_fonts_dir() -> str:
    here = os.path.dirname(__file__)
    return os.path.abspath(os.path.join(here, "..", "assets", "fonts"))

def _resolve_font_path(font_dir: str | None, font_file: str | None) -> str | None:
    if font_file:
        if os.path.isabs(font_file) and os.path.exists(font_file):
            return font_file
        if font_dir:
            p = os.path.join(font_dir, font_file)
            if os.path.exists(p):
                return p

    default_project_font = os.path.join(_project_fonts_dir(), "test.ttf")
    if os.path.exists(default_project_font):
        return default_project_font

    here = os.path.dirname(__file__)
    candidates = [p for p in [
        font_dir,
        os.path.join(here, "assets", "fonts"),
        _project_fonts_dir()
    ] if p]

    common_files = [
        "test.ttf",
        "Inter.ttf", "Inter-Regular.ttf",
        "NotoSans-Regular.ttf",
        "DejaVuSans.ttf"
    ]
    for root in candidates:
        if not os.path.isdir(root):
            continue
        for fn in common_files:
            p = os.path.join(root, fn)
            if os.path.exists(p):
                return p
    return None

def _register_font(font_path: str, base_size: int = 18):
    if dpg.does_item_exist(_FONT_REGISTRY_TAG):
        dpg.delete_item(_FONT_REGISTRY_TAG)
    with dpg.font_registry(tag=_FONT_REGISTRY_TAG):
        with dpg.font(font_path, base_size, tag=_APP_FONT_TAG):
            for hint_name in (
                "mvFontRangeHint_Default",
                "mvFontRangeHint_Cyrillic",
                "mvFontRangeHint_Korean",
                "mvFontRangeHint_Japanese",
                "mvFontRangeHint_Thai",
                "mvFontRangeHint_Vietnamese",
                "mvFontRangeHint_Chinese_Full",
            ):
                if hasattr(dpg, hint_name):
                    dpg.add_font_range_hint(getattr(dpg, hint_name))
    dpg.bind_font(_APP_FONT_TAG)
    print(f"[THEME] Loaded font: {font_path} @ {base_size}px")

# ---------------- App chrome & global ----------------
def setup_app_chrome(
    font_scale: float | None = None,
    font_dir: str | None = None,
    base_size: int = 18,
    font_file: str | None = None,
):
    dpg.configure_app(docking=True, docking_space=True)
    dpg.set_global_font_scale(font_scale if font_scale else FONT_SCALE)

    if dpg.does_item_exist(_APP_FONT_TAG):
        try:
            dpg.bind_font(_APP_FONT_TAG)
            return
        except Exception:
            pass

    if font_dir is None and font_file is None:
        font_dir = _project_fonts_dir()
        font_file = "test.ttf"

    path = _resolve_font_path(font_dir=font_dir, font_file=font_file)
    if path:
        try:
            _register_font(path, base_size=base_size)
        except Exception as e:
            print(f"[THEME] Font load failed: {e}")
    else:
        print("[THEME] No custom font found. Using Dear ImGui default.")

def set_font_scale(scale: float):
    try:
        dpg.set_global_font_scale(float(scale))
    except Exception:
        pass

# ---------------- Global themes (static) ----------------
def build_dark_theme():
    if dpg.does_item_exist("theme_dark"):
        return
    with dpg.theme(tag="theme_dark"):
        with dpg.theme_component(dpg.mvAll):
            dpg.add_theme_style(dpg.mvStyleVar_WindowPadding, SPACING, SPACING)
            dpg.add_theme_style(dpg.mvStyleVar_ItemSpacing, SPACING, SPACING)
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 9, 7)
            dpg.add_theme_style(dpg.mvStyleVar_WindowRounding, RADIUS)
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, RADIUS)
            dpg.add_theme_style(dpg.mvStyleVar_ChildRounding, RADIUS)
            dpg.add_theme_style(dpg.mvStyleVar_PopupRounding, RADIUS)
            dpg.add_theme_style(dpg.mvStyleVar_GrabRounding, RADIUS)
            dpg.add_theme_style(dpg.mvStyleVar_TabRounding, RADIUS)

            dpg.add_theme_color(dpg.mvThemeCol_Text, COLOR_TEXT_DARK)
            dpg.add_theme_color(dpg.mvThemeCol_TextDisabled, (160, 165, 170))
            dpg.add_theme_color(dpg.mvThemeCol_WindowBg, COLOR_BG_DARK_WINDOW)
            dpg.add_theme_color(dpg.mvThemeCol_ChildBg, COLOR_BG_DARK_CHILD)
            dpg.add_theme_color(dpg.mvThemeCol_PopupBg, COLOR_BG_DARK_CHILD)
            dpg.add_theme_color(dpg.mvThemeCol_FrameBg, COLOR_BG_DARK_FRAME)
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgHovered, (58, 62, 64))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgActive, (66, 70, 72))
            dpg.add_theme_color(dpg.mvThemeCol_Border, (70, 74, 78))
            dpg.add_theme_color(dpg.mvThemeCol_Separator, (70, 74, 78))

            dpg.add_theme_color(dpg.mvThemeCol_TitleBg, (18, 48, 42))
            dpg.add_theme_color(dpg.mvThemeCol_TitleBgActive, (24, 72, 63))
            dpg.add_theme_color(dpg.mvThemeCol_Tab, (28, 52, 47))
            dpg.add_theme_color(dpg.mvThemeCol_TabHovered, (34, 66, 60))
            dpg.add_theme_color(dpg.mvThemeCol_TabActive, (26, 58, 52))

def build_light_theme():
    if dpg.does_item_exist("theme_light"):
        return
    with dpg.theme(tag="theme_light"):
        with dpg.theme_component(dpg.mvAll):
            dpg.add_theme_style(dpg.mvStyleVar_WindowPadding, SPACING, SPACING)
            dpg.add_theme_style(dpg.mvStyleVar_ItemSpacing, SPACING, SPACING)
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 9, 7)
            dpg.add_theme_style(dpg.mvStyleVar_WindowRounding, RADIUS)
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, RADIUS)
            dpg.add_theme_style(dpg.mvStyleVar_ChildRounding, RADIUS)
            dpg.add_theme_style(dpg.mvStyleVar_PopupRounding, RADIUS)
            dpg.add_theme_style(dpg.mvStyleVar_GrabRounding, RADIUS)
            dpg.add_theme_style(dpg.mvStyleVar_TabRounding, RADIUS)

            dpg.add_theme_color(dpg.mvThemeCol_Text, COLOR_TEXT_LIGHT)
            dpg.add_theme_color(dpg.mvThemeCol_TextDisabled, (110, 120, 128))
            dpg.add_theme_color(dpg.mvThemeCol_WindowBg, COLOR_BG_LIGHT_WINDOW)
            dpg.add_theme_color(dpg.mvThemeCol_ChildBg, COLOR_BG_LIGHT_CHILD)
            dpg.add_theme_color(dpg.mvThemeCol_PopupBg, COLOR_BG_LIGHT_CHILD)
            dpg.add_theme_color(dpg.mvThemeCol_FrameBg, COLOR_BG_LIGHT_FRAME)
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgHovered, (218, 222, 226))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgActive, (210, 214, 219))
            dpg.add_theme_color(dpg.mvThemeCol_Border, (205, 210, 215))
            dpg.add_theme_color(dpg.mvThemeCol_Separator, (205, 210, 215))

            dpg.add_theme_color(dpg.mvThemeCol_TitleBg, (210, 232, 224))
            dpg.add_theme_color(dpg.mvThemeCol_TitleBgActive, (194, 220, 212))
            dpg.add_theme_color(dpg.mvThemeCol_Tab, (220, 234, 228))
            dpg.add_theme_color(dpg.mvThemeCol_TabHovered, (208, 226, 220))
            dpg.add_theme_color(dpg.mvThemeCol_TabActive, (200, 220, 214))

def apply_theme(dark: bool = True):
    dpg.bind_theme("theme_dark" if dark else "theme_light")

# ---------------- Component themes (static baseline) ----------------
def build_component_themes():
    if not dpg.does_item_exist("theme_button_primary"):
        with dpg.theme(tag="theme_button_primary"):
            with dpg.theme_component(dpg.mvButton):
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, RADIUS)
                dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 12, 8)
                dpg.add_theme_color(dpg.mvThemeCol_Button,        (*COLOR_PRIMARY, 255))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (*COLOR_PRIMARY_H, 255))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonActive,  (*COLOR_PRIMARY_A, 255))
                dpg.add_theme_color(dpg.mvThemeCol_Text,          (255, 255, 255))

    if not dpg.does_item_exist("theme_button_secondary"):
        with dpg.theme(tag="theme_button_secondary"):
            with dpg.theme_component(dpg.mvButton):
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, RADIUS)
                dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 10, 6)
                dpg.add_theme_color(dpg.mvThemeCol_Button, (72, 78, 86))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (88, 95, 104))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (98, 106, 116))
                dpg.add_theme_color(dpg.mvThemeCol_Text, COLOR_TEXT_DARK)

    if not dpg.does_item_exist("theme_input"):
        with dpg.theme(tag="theme_input"):
            for comp in (dpg.mvInputText, dpg.mvSliderInt, dpg.mvSliderFloat, dpg.mvCheckbox, dpg.mvCombo):
                with dpg.theme_component(comp):
                    dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, RADIUS)
                    dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 8, 6)
                    dpg.add_theme_color(dpg.mvThemeCol_FrameBg, COLOR_BG_DARK_FRAME)
                    dpg.add_theme_color(dpg.mvThemeCol_FrameBgHovered, (58, 62, 64))
                    dpg.add_theme_color(dpg.mvThemeCol_FrameBgActive, (66, 70, 72))
                    dpg.add_theme_color(dpg.mvThemeCol_Text, COLOR_TEXT_DARK)

    if not dpg.does_item_exist("theme_dialog"):
        with dpg.theme(tag="theme_dialog"):
            with dpg.theme_component(dpg.mvFileDialog):
                dpg.add_theme_style(dpg.mvStyleVar_WindowPadding, SPACING, SPACING)
                dpg.add_theme_style(dpg.mvStyleVar_WindowRounding, RADIUS)
                dpg.add_theme_color(dpg.mvThemeCol_WindowBg, (20, 21, 24))
                dpg.add_theme_color(dpg.mvThemeCol_TitleBg, (18, 48, 42))
                dpg.add_theme_color(dpg.mvThemeCol_TitleBgActive, (24, 72, 63))

# ---------------- Layout helpers ----------------
def add_spacer(height: int = 10):
    dpg.add_spacer(height=height)

@contextmanager
def hgroup(spacing: int = SPACING):
    with dpg.group(horizontal=True, horizontal_spacing=spacing):
        yield

@contextmanager
def vgroup():
    with dpg.group(horizontal=False):
        yield

@contextmanager
def card(title: str, width: int = CARD_W, height: int = CARD_H, title_color=COLOR_TITLE):
    with dpg.child_window(width=width, height=height, border=True):
        dpg.add_text(title, color=title_color)
        add_spacer(6)
        yield

def tooltip(parent_item, text: str, color=None):
    if parent_item and dpg.does_item_exist(parent_item):
        with dpg.tooltip(parent_item):
            dpg.add_text(text) if not color else dpg.add_text(text, color=color)

def bullet(text: str, parent=None):
    msg = f"• {text}"
    return dpg.add_text(msg) if parent is None else dpg.add_text(msg, parent=parent)

def header_palette():
    return {
        "fg": (210, 230, 210),
        "accent": COLOR_NOTE,
        "muted": COLOR_MUTED,
    }

# ===================== Runtime accent themes =========================
# ملاحظة: هنستخدم tag فريد في كل مرة لتجنب "Alias already exists"
_RUNTIME_COUNTERS = {"primary": 0, "sidebar": 0}
_PRIMARY_THEME_CURRENT = "theme_button_primary"   # fallback to static baseline
_SIDEBAR_ACTIVE_THEME_CURRENT = None             # اختياري (مش مستخدم حالياً)

def get_primary_button_theme_tag() -> str:
    """Used by main_ui to bind active buttons (nav/current)."""
    return _PRIMARY_THEME_CURRENT or "theme_button_primary"

def _unique_tag(prefix: str) -> str:
    _RUNTIME_COUNTERS[prefix] += 1
    return f"__rt_{prefix}_btn_{_RUNTIME_COUNTERS[prefix]}"

def _mix(c1, c2, t):
    return int(c1 + (c2 - c1) * t)

def _lighten(rgb, amount=0.25):
    return tuple(_mix(c, 255, amount) for c in rgb)

def _darken(rgb, amount=0.35):
    return tuple(_mix(c, 0, amount) for c in rgb)

def _safe_delete(tag: str):
    try:
        if tag and dpg.does_item_exist(tag):
            dpg.delete_item(tag)
    except Exception:
        # لو حصل خطأ أثناء الحذف، نتجاهله ونطلع tag جديد لاحقاً
        pass

def set_accent_color(base_rgb, hover_rgb=None, active_rgb=None):
    """
    Rebuild PRIMARY runtime theme with a fresh unique tag.
    This avoids reusing an existing tag => no alias collision.
    """
    global _PRIMARY_THEME_CURRENT

    base = tuple(base_rgb)
    hov  = tuple(hover_rgb)  if hover_rgb  else _lighten(base, 0.25)
    act  = tuple(active_rgb) if active_rgb else _darken(base, 0.35)

    # حاول احذف الـ tag القديم (لو كان runtime)
    old_tag = _PRIMARY_THEME_CURRENT if _PRIMARY_THEME_CURRENT and _PRIMARY_THEME_CURRENT.startswith("__rt_") else None
    _safe_delete(old_tag)

    # ابنِ theme جديد بوسم فريد
    new_tag = _unique_tag("primary")
    with dpg.theme(tag=new_tag):
        with dpg.theme_component(dpg.mvButton):
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, RADIUS)
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 12, 8)
            dpg.add_theme_color(dpg.mvThemeCol_Button,        (*base, 255))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (*hov, 255))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive,  (*act, 255))
            dpg.add_theme_color(dpg.mvThemeCol_Text,          (255, 255, 255))

    _PRIMARY_THEME_CURRENT = new_tag
    print(f"[THEME] Primary accent updated -> {new_tag} rgb={base}/{hov}/{act}")
