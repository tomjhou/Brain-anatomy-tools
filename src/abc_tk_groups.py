"""The "Choose groups" window of the project browser: the taxonomy as a tree (class > subclass > supertype > cluster) with a check box at every
level, for Specific group(s). Checking a group makes it one painted unit: all its cells get one color. Checking its subgroups instead (or pressing
Split) gives each of them a color of its own. Where units overlap, the cells take the finer one, so a subclass can be one color with one of its
supertypes carved out in another.

The window changes the list of units the browser keeps as text in its Group number(s) box (see UmapWindow.parse_units and set_units): the tree and the
box are two views of the same selection."""
import tkinter as tk
from tkinter import ttk, colorchooser

import numpy as np

HIERARCHY = ("class", "subclass", "supertype", "cluster")
BOX = 14                                   # the check box image's size (pixels)
MAX_ROWS = 2500                            # a search or "checked only" listing stops here
SEARCH_MS = 250


def _mix(rgb, fraction, base=(255, 255, 255)):
    return tuple(int(round(base[i] + (rgb[i] - base[i]) * fraction)) for i in range(3))


def _hex(rgb):
    return "#%02x%02x%02x" % tuple(int(v) for v in rgb)


class GroupTree:
    def __init__(self, owner):
        self.owner = owner
        self.levels = [k for k in HIERARCHY if k in owner.levels]
        self.images = {}
        self.build_taxonomy()
        top = self.top = tk.Toplevel(owner.root)
        top.title("Choose groups")
        top.geometry("640x720")
        top.minsize(420, 320)
        top.protocol("WM_DELETE_WINDOW", self.close)
        self.place()
        bar = ttk.Frame(top, padding=(8, 8, 8, 2))
        bar.pack(fill=tk.X)
        ttk.Label(bar, text="Search:").pack(side=tk.LEFT)
        self.var_search = tk.StringVar(top, "")
        entry = ttk.Entry(bar, textvariable=self.var_search)
        entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=6)
        self.var_checked = tk.BooleanVar(top, False)
        ttk.Checkbutton(bar, text="Checked only", variable=self.var_checked, command=self.rebuild).pack(side=tk.LEFT)
        self.search_job = None
        self.var_search.trace_add("write", lambda *a: self.search_soon())
        ttk.Label(top, padding=(10, 0, 10, 4), foreground="gray30", wraplength=600, justify=tk.LEFT, text=(
            "Click a box to color a group (all its cells one color). Open a group with its arrow, or press Split to give each of its subgroups its own color. "
            "Where groups overlap, the cells take the finer one.")).pack(anchor="w")
        body = ttk.Frame(top, padding=(8, 0))
        body.pack(fill=tk.BOTH, expand=True)
        self.tree = ttk.Treeview(body, show="tree", selectmode="browse")
        sb = ttk.Scrollbar(body, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.tree.tag_configure("count", foreground="gray40")
        self.tree.bind("<<TreeviewOpen>>", self.on_open)
        self.tree.bind("<ButtonRelease-1>", self.on_click)
        self.tree.bind("<space>", lambda e: (self.toggle_selected(), "break")[1])
        buttons = ttk.Frame(top, padding=8)
        buttons.pack(fill=tk.X)
        for text, command, tip in (
                ("Split", self.split_selected, "Replace the selected colored group by its subgroups, each with its own color."),
                ("Merge", self.merge_selected, "Color the selected group as one, in place of whatever is colored below it."),
                ("Color...", self.color_selected, "Choose the color of the selected colored group."),
                ("Clear all", self.clear_all, "Remove every group."),
                ("Close", self.close, "Close this window (the same as its X).")):
            b = ttk.Button(buttons, text=text, command=command)
            b.pack(side=tk.LEFT, padx=(0, 6))
        self.status = ttk.Label(buttons, foreground="gray30")
        self.status.pack(side=tk.LEFT, padx=8)
        self.rebuild(expand_painted=True)
        entry.focus_set()

    # ---- the taxonomy: which groups there are, and what is under what

    def build_taxonomy(self):
        owner = self.owner
        self.nodes = {}                 # (level, key) -> {"label", "count", "level", "key"}
        self.parent = {}                # (level, key) -> (level, key) of the group above
        self.children = {}              # (level, key) -> [the groups below]
        per_level = []
        for level in self.levels:
            codes, uniques = owner.category_codes(level)
            keys = owner.group_keys(level, uniques)
            counts = np.bincount(codes[codes >= 0], minlength=len(uniques))
            keyed = []
            for u in range(len(uniques)):
                if counts[u] == 0:
                    continue
                try:
                    key = int(keys[u])
                except (TypeError, ValueError):
                    continue
                self.nodes[(level, key)] = {"label": str(uniques[u]), "count": int(counts[u]), "level": level, "key": key}
                keyed.append(key)
            per_level.append((level, codes, keys, len(uniques)))
        for (upper, codes_up, keys_up, n_up), (lower, codes_low, keys_low, n_low) in zip(per_level, per_level[1:]):
            both = (codes_up >= 0) & (codes_low >= 0)
            pairs = np.unique(codes_up[both].astype(np.int64) * n_low + codes_low[both])
            for p in pairs:
                up, low = int(p // n_low), int(p % n_low)
                try:
                    a, b = (upper, int(keys_up[up])), (lower, int(keys_low[low]))
                except (TypeError, ValueError):
                    continue
                if a in self.nodes and b in self.nodes and b not in self.parent:
                    self.parent[b] = a
                    self.children.setdefault(a, []).append(b)
        for kids in self.children.values():
            kids.sort(key=lambda n: n[1])
        top_level = self.levels[0] if self.levels else None
        self.roots = sorted((n for n in self.nodes if n[0] == top_level), key=lambda n: n[1]) if top_level else []

    def ancestors(self, node):
        while node in self.parent:
            node = self.parent[node]
            yield node

    # ---- what is painted

    def units(self):
        return self.owner.parse_units()[0] if self.owner.mode == "single" else []

    def states(self):
        """(painted: {node: rgb}, covered: {node: rgb of the painted group above}, partial: set of nodes with something painted below)"""
        units = self.units()
        painted, partial = {}, set()
        for unit in units:
            if unit in self.nodes:
                painted[unit] = tuple(int(v) for v in self.owner.unit_rgb(unit, units))
        covered = {}
        for node in self.nodes:
            for above in self.ancestors(node):
                if above in painted:
                    covered[node] = painted[above]
                    break
        for node in painted:
            partial.update(self.ancestors(node))
        return painted, covered, partial

    def box_image(self, kind, rgb):
        rgb = tuple(int(v) for v in rgb) if rgb is not None else None
        key = (kind, rgb)
        img = self.images.get(key)
        if img is not None:
            return img
        img = tk.PhotoImage(master=self.top, width=BOX, height=BOX)
        border = "#6b6b6b" if kind in ("on", "off") else "#b0b0b0"
        if kind == "on":
            fill = _hex(rgb)
        elif kind == "covered":
            fill = _hex(_mix(rgb, 0.35))
        else:
            fill = "#ffffff"
        for y in range(BOX):
            row = []
            for x in range(BOX):
                edge = x in (0, BOX - 1) or y in (0, BOX - 1)
                color = border if edge else fill
                if kind == "partial" and not edge and 4 <= x <= BOX - 5 and 4 <= y <= BOX - 5:
                    color = _hex(rgb)
                row.append(color)
            img.put("{" + " ".join(row) + "}", to=(0, y))
        if kind == "on":
            ink = "#ffffff" if sum(rgb) < 460 else "#202020"
            for x, y in ((3, 7), (4, 8), (5, 9), (6, 8), (7, 7), (8, 6), (9, 5), (10, 4), (3, 6), (4, 7), (5, 8), (6, 7), (7, 6), (8, 5), (9, 4)):
                img.put(ink, to=(x, y))
        self.images[key] = img
        return img

    def node_text(self, node):
        info = self.nodes[node]
        return f"{info['label']}   ({info['count']:,})"

    def row_image(self, node, painted, covered, partial):
        if node in painted:
            return self.box_image("on", painted[node])
        if node in covered:
            return self.box_image("covered", covered[node])
        if node in partial:
            units = self.units()
            below = next((u for u in units if u in self.nodes and node in set(self.ancestors(u))), None)
            return self.box_image("partial", tuple(int(v) for v in self.owner.unit_rgb(below, units)) if below is not None else (120, 120, 120))
        return self.box_image("off", None)

    # ---- the tree on the screen

    def iid(self, node):
        return f"{node[0]}|{node[1]}"

    def node_of(self, iid):
        level, key = iid.split("|")
        return (level, int(key))

    def insert(self, parent_iid, node, children_below, painted, covered, partial, opened=False):
        iid = self.iid(node)
        if self.tree.exists(iid):
            return iid
        self.tree.insert(parent_iid, "end", iid=iid, text=self.node_text(node), image=self.row_image(node, painted, covered, partial), open=opened)
        if children_below and self.children.get(node):
            self.tree.insert(iid, "end", iid=iid + "|...", text="")          # (a placeholder, replaced by the children when it is opened)
        return iid

    def rebuild(self, expand_painted=False):
        painted, covered, partial = self.states()
        self.tree.delete(*self.tree.get_children())
        query = self.var_search.get().strip().casefold()
        if query or self.var_checked.get():
            shown = set()
            for node, info in self.nodes.items():
                if self.var_checked.get():
                    if node not in painted:
                        continue
                if query and query not in info["label"].casefold():
                    continue
                shown.add(node)
                shown.update(self.ancestors(node))
            order = sorted(shown, key=lambda n: (self.levels.index(n[0]), n[1]))
            capped = len(order) > MAX_ROWS
            for node in order[:MAX_ROWS]:
                parent = self.parent.get(node)
                self.insert(self.iid(parent) if parent in shown and self.tree.exists(self.iid(parent)) else "", node, False, painted, covered, partial, opened=True)
            self.status.config(text=f"{len(order):,} groups" + (f" (the first {MAX_ROWS:,} shown)" if capped else ""))
            return
        for node in self.roots:
            self.insert("", node, True, painted, covered, partial)
        self.status.config(text="")
        if expand_painted:
            for unit in painted:
                self.reveal(unit)

    def reveal(self, node):
        chain = list(self.ancestors(node))[::-1]
        for above in chain:
            iid = self.iid(above)
            if self.tree.exists(iid):
                self.open_node(iid)
        iid = self.iid(node)
        if self.tree.exists(iid):
            self.tree.see(iid)

    def open_node(self, iid):
        """Put the real children under a row that has only its placeholder, and open it."""
        placeholder = iid + "|..."
        if self.tree.exists(placeholder):
            self.tree.delete(placeholder)
            painted, covered, partial = self.states()
            for kid in self.children.get(self.node_of(iid), []):
                self.insert(iid, kid, True, painted, covered, partial)
        self.tree.item(iid, open=True)

    def on_open(self, event):
        iid = self.tree.focus()
        if iid and not iid.endswith("|..."):
            self.open_node(iid)

    def refresh(self):
        """The check boxes and colors again (the selection changed: a legend click, the box, Load state, or this window)."""
        painted, covered, partial = self.states()
        stack = list(self.tree.get_children())
        while stack:
            iid = stack.pop()
            stack.extend(self.tree.get_children(iid))
            if iid.endswith("|..."):
                continue
            self.tree.item(iid, image=self.row_image(self.node_of(iid), painted, covered, partial))

    # ---- clicks

    def on_click(self, event):
        if self.tree.identify_region(event.x, event.y) not in ("tree", "cell"):
            return
        if self.tree.identify_element(event.x, event.y) != "image":
            return
        iid = self.tree.identify_row(event.y)
        if iid and not iid.endswith("|..."):
            self.toggle(self.node_of(iid))

    def selected(self):
        iid = self.tree.focus() or (self.tree.selection() or [""])[0]
        return self.node_of(iid) if iid and not iid.endswith("|...") else None

    def toggle_selected(self):
        node = self.selected()
        if node is not None:
            self.toggle(node)

    def toggle(self, node):
        units = list(self.units())
        if node in units:
            units.remove(node)
        else:
            units.append(node)
        self.apply(units)

    def split_selected(self):
        node = self.selected()
        units = list(self.units())
        kids = self.children.get(node, []) if node is not None else []
        if node is None or not kids:
            self.status.config(text="Select a group that has subgroups.")
            return
        if node not in units:
            units.append(node)                          # (splitting a group that is not colored yet: its subgroups, each in a color)
        at = units.index(node)
        units.remove(node)
        for kid in kids:
            if kid not in units:
                units.insert(at, kid)
                at += 1
        self.apply(units)
        self.open_node(self.iid(node))

    def merge_selected(self):
        node = self.selected()
        if node is None:
            self.status.config(text="Select the group to color as one.")
            return
        units = [u for u in self.units() if not (u in self.nodes and node in set(self.ancestors(u))) and u != node]
        units.append(node)
        self.apply(units)

    def color_selected(self):
        node = self.selected()
        units = self.units()
        if node is None or node not in units:
            self.status.config(text="Select a colored group (one with a filled box) to change its color.")
            return
        current = self.owner.unit_rgb(node, units)
        chosen = colorchooser.askcolor(color=_hex(current), parent=self.top, title="Color of " + self.nodes[node]["label"])
        if chosen and chosen[0]:
            self.owner.set_unit_color(node, tuple(int(round(v)) for v in chosen[0]))

    def clear_all(self):
        self.apply([u for u in self.units() if u not in self.nodes])

    def apply(self, units):
        self.status.config(text="")
        self.owner.set_units(units)

    # ---- search

    def search_soon(self):
        if self.search_job is not None:
            self.top.after_cancel(self.search_job)
        self.search_job = self.top.after(SEARCH_MS, self.rebuild)

    # ---- the window

    def place(self):
        owner = self.owner
        try:
            owner.root.update_idletasks()
            x = owner.root.winfo_rootx() + max(40, owner.root.winfo_width() - 700)
            y = owner.root.winfo_rooty() + 60
            self.top.geometry(f"+{x}+{y}")
        except tk.TclError:
            pass

    def lift(self):
        try:
            self.top.deiconify()
            self.top.lift()
            self.top.focus_force()
        except tk.TclError:
            pass

    def close(self):
        self.owner.groups_window = None
        try:
            self.top.destroy()
        except tk.TclError:
            pass
