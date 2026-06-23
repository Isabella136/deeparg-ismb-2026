import math

import numpy as np
from matplotlib import axes
from matplotlib import text as mtext
from matplotlib.backend_bases import RendererBase

# Thanking the gods at stackoverflow for this:
# https://stackoverflow.com/questions/19353576/curved-text-rendering-in-matplotlib


class CurvedText(mtext.Text):
    """A text object that follows an arbitrary curve."""

    x: np.array
    y: np.array
    zorder: float
    characters: tuple[(str, mtext.Text)]
    va: str

    def __init__(
        self, x: np.array, y: np.array, text: str, ax: axes.Axes, **kwargs  # noqa: ANN003
    ) -> None:
        """Initialize curved text."""
        super().__init__(x[0], y[0], " ", **kwargs)

        ax.add_artist(self)

        # saving the curve:
        self.x = x
        self.y = y
        self.zorder = self.get_zorder()
        self.va = self.get_verticalalignment()

        # creating the text objects
        self.characters = []
        for c in text:
            if c == " ":
                # make this an invisible 'a':
                t = mtext.Text(0, 0, "a")
                t.set_alpha(0.0)
            else:
                t = mtext.Text(0, 0, c, **kwargs)

            # resetting unnecessary arguments
            t.set_ha("center")
            t.set_rotation(0)
            t.set_zorder(self.zorder + 1)

            self.characters.append((c, t))
            ax.add_artist(t)

    def set_zorder(self, zorder: float) -> None:
        """Set the zorder for the artist.  Artists with lower zorder values are drawn first."""  # noqa: E501
        super().set_zorder(zorder)
        self.zorder = self.get_zorder()
        for c, t in self.characters:
            t.set_zorder(self.zorder + 1)

    def draw(self, renderer: RendererBase, *args, **kwargs) -> None:  # noqa: ANN002, ANN003, ARG002
        """Overload Text.draw() to update the positions and rotation angles of self.characters."""  # noqa: E501
        self.make_text_centered(renderer)
        self.make_text_oriented()
        self.update_positions(renderer)

    def get_aspect_ratio(self) -> float:
        """Get the figure's aspect ratio."""
        # From https://stackoverflow.com/questions/41597177/get-aspect-ratio-of-axes/42014041#42014041

        # Get x and y limits
        xlim = self.axes.get_xlim()
        ylim = self.axes.get_ylim()

        # Get total figure size in inches
        fig_w, fig_h = self.axes.get_figure().get_size_inches()

        # Get axis size on figure
        _, _, ax_w, ax_h = self.axes.get_position().bounds

        # Find ratio of display units
        disp_ratio = (fig_h * ax_h) / (fig_w * ax_w)

        # Find ratio of limit units
        limit_ratio = (ylim[1] - ylim[0]) / (xlim[1] - xlim[0])

        # Return 1 over aspect, which is disp_ratio / limit_ratio
        return 1 / (disp_ratio / limit_ratio)

    def get_arc_data(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Retrieve arc data information."""
        # Find points of the curve in figure coordinates
        x_fig, y_fig = (
            np.array(l)
            for l in zip(  # noqa: E741
                *self.axes.transData.transform([
                    (i, j) for i, j in zip(self.x, self.y)
                ])
            )
        )

        # Calculate point distances
        x_fig_dist = x_fig[1:] - x_fig[:-1]
        y_fig_dist = y_fig[1:] - y_fig[:-1]

        r_fig_dist = np.sqrt(x_fig_dist**2 + y_fig_dist**2)

        # Calculate cumulative arc length in figure coordinates
        l_fig = np.insert(np.cumsum(r_fig_dist), 0, 0)

        # Find angles in figure coordinates
        rads = np.arctan2((y_fig[1:] - y_fig[:-1]), (x_fig[1:] - x_fig[:-1]))
        degs = np.rad2deg(rads)

        return (r_fig_dist, l_fig, rads, degs)

    def make_text_centered(self, renderer: RendererBase) -> None:
        """Curtail x and y array to make text centered."""
        # Find the arc length
        _, l_fig, _, _ = self.get_arc_data()
        arc_length = l_fig[-1]

        # This is total length of text
        text_length = 0

        for c, t in self.characters:
            # Figure out character width as well as optimal start
            t.set_rotation(0)
            t.set_va("center")
            bbox1 = t.get_window_extent(renderer=renderer)
            w = bbox1.width

            # finding the two data points between which the horizontal
            # center point of the character will be situated
            # left and right indices:
            il = np.where(text_length + w / 2 >= l_fig)[0][-1]

            # how much of the letter width was needed to find il:
            used = l_fig[il] - text_length
            text_length = l_fig[il]

            # updating rel_pos to right edge of character
            text_length += w - used - 2

        # Get the start and end of text
        text_start = arc_length / 2 - text_length / 2
        text_end = arc_length / 2 + text_length / 2

        # Figure out the left- and right-most indices
        il = np.where(text_start <= l_fig)[0][0]
        ir = np.where(text_end >= l_fig)[0][-1]

        # Update x and y
        self.x = self.x[il:ir + 1]
        self.y = self.y[il:ir + 1]

    def make_text_oriented(self) -> None:
        """Will reverse character order if text is in bottom half of figure."""
        y_mid = np.median(self.y)
        if y_mid < 0:
            self.x = np.flip(self.x)
            self.y = np.flip(self.y)
            if self.va == "bottom":
                self.va = "top"
            elif self.va == "top":
                self.va = "bottom"

    def update_positions(self, renderer: RendererBase) -> None:  # noqa: PLR0914
        """Update positions and rotations of the individual text elements."""
        # Determine the aspect ratio
        aspect_ratio = self.get_aspect_ratio()

        # Get arc distance, cummulative length, and angles in radian and degree
        r_fig_dist, l_fig, rads, degs = self.get_arc_data()

        rel_pos = 0
        for c, t in self.characters:
            # finding the width of c:
            t.set_rotation(0)
            t.set_va("center")
            bbox1 = t.get_window_extent(renderer=renderer)
            w = bbox1.width

            # ignore all letters that don't fit:
            if rel_pos + w / 2 > l_fig[-1]:
                t.set_alpha(0.0)
                rel_pos += w
                continue

            if c != " ":
                t.set_alpha(1.0)

            # finding the two data points between which the horizontal
            # center point of the character will be situated
            # left and right indices:
            il = np.where(rel_pos + w / 2 >= l_fig)[0][-1]
            ir = np.where(rel_pos + w / 2 <= l_fig)[0][0]

            # if we exactly hit a data point:
            if ir == il:
                ir += 1

            # how much of the letter width was needed to find il:
            used = l_fig[il] - rel_pos
            rel_pos = l_fig[il]

            # relative distance between il and ir where the center
            # of the character will be
            fraction = (w / 2 - used) / r_fig_dist[il]

            # setting the character position in data coordinates:
            # interpolate between the two points:
            x = self.x[il] + fraction * (self.x[ir] - self.x[il])
            y = self.y[il] + fraction * (self.y[ir] - self.y[il])

            # getting the offset when setting correct vertical alignment
            # in data coordinates
            t.set_va(self.va)
            bbox2 = t.get_window_extent(renderer=renderer)

            bbox1d = self.axes.transData.inverted().transform(bbox1)
            bbox2d = self.axes.transData.inverted().transform(bbox2)
            dr = np.array(bbox2d[0] - bbox1d[0])

            # the rotation/stretch matrix
            rad = rads[il]
            rot_mat = np.array([
                [math.cos(rad), math.sin(rad) * aspect_ratio],
                [-math.sin(rad) / aspect_ratio, math.cos(rad)],
            ])

            # computing the offset vector of the rotated character
            drp = np.dot(dr, rot_mat)

            # setting final position and rotation:
            t.set_position(np.array([x, y]) + drp)
            t.set_rotation(degs[il])

            t.set_va("center")
            t.set_ha("center")

            # updating rel_pos to right edge of character
            rel_pos += w - used - 2
