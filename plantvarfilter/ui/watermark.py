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

    if hasattr(dpg, "get_item_rect_min") and hasattr(dpg, "get_item_rect_max"):
        try:
            return dpg.get_item_rect_min(tag), dpg.get_item_rect_max(tag)
        except Exception:
            pass

    try:
        pos = dpg.get_item_pos(tag)
        w, h = dpg.get_item_rect_size(tag)
        return pos, (pos[0] + w, pos[1] + h)
    except Exception:
        pass

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
    """Draw the main watermark (assets/logo.png) centered inside target window."""
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
                    width: int = 220,
                    margin=(20, 20)):
    """Place the lab logo at bottom-right of the workspace with right-aligned text on its left.
       Repositions existing items if they already exist (responsive)."""
    here = os.path.dirname(__file__)
    img_path = os.path.abspath(os.path.join(here, "..", "assets", image_name))

    tex_wh = _ensure_texture(_LAB_TEX, img_path)
    if not tex_wh or not dpg.does_item_exist(target_window_tag):
        return
    tex_w, tex_h = tex_wh

    # find top-level window to attach items to
    def _top_window(tag: str):
        cur = tag
        for _ in range(20):
            info = dpg.get_item_info(cur)
            if not info:
                break
            if "Window" in info.get("type", ""):
                return cur
            parent = info.get("parent")
            if not parent:
                break
            cur = parent
        return None

    host_window = _top_window(target_window_tag)
    if not host_window:
        return

    # absolute rects
    pmin_ca, pmax_ca = _window_abs_rect(target_window_tag)
    pmin_win, _ = _window_abs_rect(host_window)

    # sizes
    ca_w = int(pmax_ca[0] - pmin_ca[0])
    target_w = int(min(width if width else ca_w * 0.18, tex_w))
    target_h = max(1, int(target_w * (tex_h / max(tex_w, 1))))

    # logo pos (abs -> local)
    abs_logo_x = int(pmax_ca[0] - target_w - int(margin[0]))
    abs_logo_y = int(pmax_ca[1] - target_h - int(margin[1]))
    logo_pos = (abs_logo_x - int(pmin_win[0]), abs_logo_y - int(pmin_win[1]))

    # text (right-aligned block to the left)
    main_text = "Ye-Lab, PKU-IAAS"
    sub_text  = ("Develop by:"
                 " Ahmed Yassin and  Falak Sher Khan")
    gap = 12
    main_w, main_h = dpg.get_text_size(main_text)
    sub_w,  sub_h  = dpg.get_text_size(sub_text)
    text_block_w = max(main_w, sub_w)
    text_block_h = main_h + 4 + sub_h

    text_x = logo_pos[0] - gap - text_block_w
    text_y = int(logo_pos[1] + (target_h - text_block_h) / 2)

    # clamp to keep text inside the content area (left edge)
    ca_local_x0 = int(pmin_ca[0] - pmin_win[0])
    text_x = max(ca_local_x0 + 8, text_x)

    # ensure a container group (so we can move everything together)
    if not dpg.does_item_exist(_LAB_IMG):
        dpg.add_group(parent=host_window, tag=_LAB_IMG)

    # logo item
    _ICON = "__lab_icon__"
    if dpg.does_item_exist(_ICON):
        dpg.configure_item(_ICON, pos=logo_pos, width=target_w, height=target_h, parent=_LAB_IMG)
    else:
        dpg.add_image(_LAB_TEX, pos=logo_pos, width=target_w, height=target_h, parent=_LAB_IMG, tag=_ICON)

    # text items (right-aligned under each other)
    _T1, _T2 = "__lab_t1__", "__lab_t2__"
    if dpg.does_item_exist(_T1):
        dpg.configure_item(_T1, pos=(text_x + (text_block_w - main_w), text_y), parent=_LAB_IMG)
        dpg.configure_item(_T2, pos=(text_x + (text_block_w - sub_w),  text_y + main_h + 4), parent=_LAB_IMG)
    else:
        dpg.add_text(main_text, pos=(text_x + (text_block_w - main_w), text_y),
                     color=(255, 255, 255, 255), parent=_LAB_IMG, tag=_T1)
        dpg.add_text(sub_text,  pos=(text_x + (text_block_w - sub_w),  text_y + main_h + 4),
                     color=(255, 255, 255, 255), parent=_LAB_IMG, tag=_T2)



    print(f"[wm] Logo at {logo_pos} size {target_w}x{target_h}; text at ({text_x},{text_y}) width {text_block_w}")


