import os
import dearpygui.dearpygui as dpg

# viewport overlay drawlist
_VP_DL = "__wm_vp_drawlist__"

# main watermark (center)
_TEX = "__wm_texture__"
_IMG = "__wm_image__"

# lab signature (bottom-right)
_LAB_TEX = "__wm_lab_texture__"
_LAB_IMG = "__wm_lab_image__"


def _ensure_viewport_drawlist(front: bool = True):
    if not dpg.does_item_exist(_VP_DL):
        dpg.add_viewport_drawlist(front=front, tag=_VP_DL)
    else:
        dpg.configure_item(_VP_DL, front=front)


def _window_abs_rect(tag: str):
    """Return absolute (pmin, pmax) for `tag` or fallback to viewport size."""
    if not dpg.does_item_exist(tag):
        return (0, 0), (dpg.get_viewport_client_width(), dpg.get_viewport_client_height())

    # Preferred modern API
    if hasattr(dpg, "get_item_rect_min") and hasattr(dpg, "get_item_rect_max"):
        try:
            return dpg.get_item_rect_min(tag), dpg.get_item_rect_max(tag)
        except Exception:
            pass

    # Fallback: use pos + rect size
    try:
        pos = dpg.get_item_pos(tag)
        w, h = dpg.get_item_rect_size(tag)
        return pos, (pos[0] + w, pos[1] + h)
    except Exception:
        pass

    # Final fallback
    return (0, 0), (dpg.get_viewport_client_width(), dpg.get_viewport_client_height())


def _ensure_texture(tag: str, logo_path: str):
    """Create a static texture if not exists; return (width, height) or None."""
    if not os.path.exists(logo_path):
        print(f"[wm] logo not found: {logo_path}")
        return None
    if not dpg.does_item_exist(tag):
        w, h, c, data = dpg.load_image(logo_path)
        with dpg.texture_registry(show=False):
            dpg.add_static_texture(width=w, height=h, default_value=data, tag=tag)
        return w, h
    info = dpg.get_item_configuration(tag)
    return info.get("width", 1), info.get("height", 1)


def _delete_item_if_exists(tag: str):
    if dpg.does_item_exist(tag):
        try:
            dpg.delete_item(tag)
        except Exception:
            pass


def setup(alpha: int = 40,
          scale: float = 0.40,
          target_window_tag: str = "content_area",
          front: bool = True):
    """
    Draw the main watermark (assets/logo.png) centered inside target window.
    """
    here = os.path.dirname(__file__)
    logo_path = os.path.abspath(os.path.join(here, "..", "assets", "logo.png"))

    tex_wh = _ensure_texture(_TEX, logo_path)
    if not tex_wh:
        return
    tex_w, tex_h = tex_wh

    _ensure_viewport_drawlist(front=front)

    pmin, pmax = _window_abs_rect(target_window_tag)
    rect_w = max(0, int(pmax[0] - pmin[0]))
    rect_h = max(0, int(pmax[1] - pmin[1]))

    target_w = max(1, int(tex_w * float(scale)))
    target_h = max(1, int(tex_h * float(scale)))

    x = int(pmin[0] + (rect_w - target_w) / 2)
    y = int(pmin[1] + (rect_h - target_h) / 2)

    a = max(0, min(255, int(255 * (float(alpha) / 100.0))))

    _delete_item_if_exists(_IMG)
    dpg.draw_image(
        _TEX,
        pmin=(x, y),
        pmax=(x + target_w, y + target_h),
        uv_min=(0, 0),
        uv_max=(1, 1),
        color=(255, 255, 255, a),
        parent=_VP_DL,
        tag=_IMG,
    )
    print(f"[wm] Watermark placed @ {x},{y} size {target_w}x{target_h} alpha={a}")


def place_signature(target_window_tag: str = "content_area",
                    image_name: str = "logo_lab.png",
                    alpha: int = 90,          # ignored for opacity; kept for API compatibility
                    width: int = 220,
                    margin=(20, 20),
                    front: bool = True,
                    add_bg_card: bool = True,
                    pad: int = 10,
                    rounding: int = 10,
                    border_alpha: int = 60):  # ignored for opacity; kept for API compatibility
    """Bottom-right lab signature with an optional opaque background card (no transparency)."""
    here = os.path.dirname(__file__)
    img_path = os.path.abspath(os.path.join(here, "..", "assets", image_name))

    tex_wh = _ensure_texture(_LAB_TEX, img_path)
    if not tex_wh:
        return
    tex_w, tex_h = tex_wh

    _ensure_viewport_drawlist(front=front)

    pmin, pmax = _window_abs_rect(target_window_tag)
    rect_w = max(0, int(pmax[0] - pmin[0]))
    rect_h = max(0, int(pmax[1] - pmin[1]))

    # Do not upscale beyond native texture width to avoid blurring
    target_w = int(min(width if width else rect_w * 0.18, tex_w))
    target_h = max(1, int(target_w * (tex_h / max(tex_w, 1))))

    x = int(pmax[0] - target_w - int(margin[0]))
    y = int(pmax[1] - target_h - int(margin[1]))

    # Force full opacity for both the image and the optional card
    a_img = 255
    a_border = 255

    _delete_item_if_exists(_LAB_IMG)

    if add_bg_card:
        # Opaque card to maximize contrast on dark UI
        bg_fill = (20, 20, 20, 255)       # solid dark fill (no transparency)
        bg_border = (255, 255, 255, a_border)
        dpg.draw_rectangle(
            pmin=(x - pad, y - pad),
            pmax=(x + target_w + pad, y + target_h + pad),
            color=bg_border,
            fill=bg_fill,
            rounding=rounding,
            parent=_VP_DL,
        )

    dpg.draw_image(
        _LAB_TEX,
        pmin=(x, y),
        pmax=(x + target_w, y + target_h),
        uv_min=(0, 0),
        uv_max=(1, 1),
        color=(255, 255, 255, a_img),     # fully opaque
        parent=_VP_DL,
        tag=_LAB_IMG,
    )

    print(f"[wm] Signature placed @ {x},{y} size {target_w}x{target_h} opaque")
