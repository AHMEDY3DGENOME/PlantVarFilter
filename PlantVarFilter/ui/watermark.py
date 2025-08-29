import os
import dearpygui.dearpygui as dpg

_VP_DL = "__wm_vp_drawlist__"   # viewport overlay drawlist
_TEX   = "__wm_texture__"
_IMG   = "__wm_image__"


def _ensure_texture(logo_path: str):
    if not os.path.exists(logo_path):
        print("[wm] logo not found:", logo_path)
        return None
    if not dpg.does_item_exist(_TEX):
        w, h, c, data = dpg.load_image(logo_path)
        with dpg.texture_registry(show=False):
            dpg.add_static_texture(width=w, height=h, default_value=data, tag=_TEX)
        return w, h
    info = dpg.get_item_configuration(_TEX)
    return info.get("width", 1), info.get("height", 1)


def _ensure_viewport_drawlist(front: bool = True):
    if not dpg.does_item_exist(_VP_DL):
        dpg.add_viewport_drawlist(front=front, tag=_VP_DL)
    else:
        dpg.configure_item(_VP_DL, front=front)


def _delete_old_image():
    if dpg.does_item_exist(_IMG):
        try:
            dpg.delete_item(_IMG)
        except Exception:
            pass


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

    # Fallback: use pos + rect size when available
    try:
        pos = dpg.get_item_pos(tag)  # may be (0,0) if not explicitly set
        w, h = dpg.get_item_rect_size(tag)
        return pos, (pos[0] + w, pos[1] + h)
    except Exception:
        pass

    # Final fallback: viewport
    return (0, 0), (dpg.get_viewport_client_width(), dpg.get_viewport_client_height())


def setup(alpha: int = 40, scale: float = 0.40,
          target_window_tag: str = "content_area",
          front: bool = True):
    """
    Draw watermark centered over the absolute rect of `target_window_tag`
    as a viewport overlay (always visible).
    """
    here = os.path.dirname(__file__)
    logo_path = os.path.abspath(os.path.join(here, "..", "assets", "logo.png"))

    tex_wh = _ensure_texture(logo_path)
    if not tex_wh:
        return
    tex_w, tex_h = tex_wh

    _ensure_viewport_drawlist(front=front)

    pmin, pmax = _window_abs_rect(target_window_tag)
    rect_w = max(0, int(pmax[0] - pmin[0]))
    rect_h = max(0, int(pmax[1] - pmin[1]))

    # target size
    target_w = max(1, int(tex_w * float(scale)))
    target_h = max(1, int(tex_h * float(scale)))

    # center inside content_area rect
    x = int(pmin[0] + (rect_w - target_w) / 2)
    y = int(pmin[1] + (rect_h - target_h) / 2)

    a = max(0, min(255, int(255 * (float(alpha) / 100.0))))

    _delete_old_image()
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
