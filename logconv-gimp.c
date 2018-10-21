/*
 * logconv.c 1.3
 *
 * Copyright (C) 1998-2018 Alessandro Baldoni
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
/*
 * This is version 1.3
 * LogConv for RGB images works on their luminance.
 */
#include <glib.h>
#include <gtk/gtk.h>
#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>
#include <math.h>

#define KD(sigma, ag, al) ((sigma * M_PI)/ (sqrt (pow (ag, 2.0) + pow(al, 2.0))))
#define KS(sigma, kd, al) ((sigma * M_PI) / (al * kd))
#define SG(sigma, ks) (sigma * sqrt (1.0 - (1.0 / pow (ks, 2.0))))
#define SL(sigma, ks, kd) (sigma / (ks * kd))

typedef struct _entry {
	gdouble Ag, Al;
} entry;

static void ZoomImage(gdouble *, gdouble *, gint, gint, gint);
static gdouble *LoGConvolution(GimpDrawable *, gdouble, gint, gint);
static gint LoGZero(GimpDrawable *, gdouble, gint, gint, gint, gint, gint);
static gdouble *MakeFilter(gdouble, gint, gint);
static gint oddsupport(gdouble);
static void query();
static void run(const char *, int, const GimpParam *, int *, GimpParam **);
static gint LoGdialog(void);
static inline gdouble SobelX(gdouble *, gint, gint, gint);
static inline gdouble SobelY(gdouble *, gint, gint, gint);
static inline gdouble RobertsX(gdouble *, gint, gint, gint);
static inline gdouble RobertsY(gdouble *, gint, gint, gint);
static inline gdouble LuminanceSRGB(gint, gint, gint);
static inline gdouble LuminancePerceptive1(gint, gint, gint);
static inline gdouble LuminancePerceptive2(gint, gint, gint);
static inline gdouble H1(gdouble, gdouble);
static inline gdouble H2(gdouble, gdouble);
static inline gdouble grad(gdouble, gdouble);

/*
 * Precomputed values of Al e Ag for a given pa
 */
static entry Table1[11] = { { 3.554140, 4.205983 }, /* pa =  0.0001 */
{ 3.402057, 4.061595 }, /* pa =  0.0003 */
{ 3.227792, 3.896178 }, /* pa =  0.0010 */
{ 3.060844, 3.737694 }, /* pa =  0.0030 */
{ 2.867757, 3.554300 }, /* pa =  0.0100 */
{ 2.712526, 3.406712 }, /* pa =  0.0250 */
{ 2.461219, 3.167270 }, /* pa =  0.1000 */
{ 2.244686, 2.960158 }, /* pa =  0.3000 */
{ 1.984301, 2.709568 }, /* pa =  1.0000 */
{ 1.718008, 2.450782 }, /* pa =  3.0000 */
{ 1.378026, 2.115151 } /* pa = 10.0000 */
};

static char* AllowablePA[11] = {
		"0.0001","0.0003", "0.0010", "0.0030", "0.0100", "0.0250", "0.1000", "0.3000", "1.0000", "3.0000", "10.0000"
};

static char* LoGType[3] = {
		"Standard", "with Roberts", "with Sobel"
};

static char* LuminanceType[3] = {
		"sRGB", "Perceived 1", "Perceived 2"
};

typedef struct _LoGDialog {
	gint run, pa, pc1, pc2, which, luminance;
	gdouble sd;
} LoGDialog;

static LoGDialog ldial = { FALSE, 0, 25, 60, 0, 1, 2.0 };

GimpPlugInInfo PLUG_IN_INFO = {
NULL, /* init_proc */
NULL, /* quit_proc */
query, /* query_proc */
run, /* run_proc */
};

MAIN ()

static void query() {
	static GimpParamDef log_args[] =
			{ { GIMP_PDB_INT32, "run_mode", "Interactive, non-interactive" },
			  { GIMP_PDB_IMAGE, "image", "Input image" },
			  {	GIMP_PDB_DRAWABLE, "drawable", "Input drawable" },
			  {	GIMP_PDB_INT32, "pa", "Selected pa (0..10)" },
			  {	GIMP_PDB_FLOAT, "sigma", "Standard deviation" },
			  {	GIMP_PDB_INT32, "type",	"0: Standard LoG, 1: LoG with Roberts, 2: LoG with Sobel" },
			  { GIMP_PDB_INT32, "pc1", "Low threshold for Roberts & Sobel"},
			  { GIMP_PDB_INT32, "pc2", "High threshold for Roberts & Sobel"},
			  { GIMP_PDB_INT32, "luminance", "0: sRGB, 1: perceived 1, 2: perceived 2"}
			};
	static int nlog_args = sizeof(log_args) / sizeof(log_args[0]);

	gimp_install_procedure("plug_in_LoG", "Apply the LoG filter", "",
			"Alessandro Baldoni", "Alessandro Baldoni", "1998-2018",
			"<Image>/Filters/Edge-Detect/LoG", "RGB,GRAY", GIMP_PLUGIN,
			nlog_args, 0, log_args, NULL);
}

static void run(const char *name, int nparams, const GimpParam * param, int *nreturn_vals,
		GimpParam ** return_vals) {
	static GimpParam values[1];
	GimpRunMode run_mode;
	GimpDrawable *draw;

	run_mode = param[0].data.d_int32;

	*nreturn_vals = 1;
	*return_vals = values;

	values[0].type = GIMP_PDB_STATUS;
	values[0].data.d_status = GIMP_PDB_CALLING_ERROR;

	draw = gimp_drawable_get(param[2].data.d_int32);
	switch (run_mode) {
	case GIMP_RUN_INTERACTIVE:
		gimp_get_data("plug_in_LoG", &ldial);
		if (LoGdialog()) {
			if (LoGZero(draw, ldial.sd, ldial.pa, ldial.which, ldial.pc1, ldial.pc2, ldial.luminance) != -1) {
				values[0].data.d_status = GIMP_PDB_SUCCESS;
				gimp_displays_flush();
				gimp_procedural_db_set_data("plug_in_LoG", &ldial, sizeof(LoGDialog));
			} else
				values[0].data.d_status = GIMP_PDB_EXECUTION_ERROR;
		}
		break;
	case GIMP_RUN_NONINTERACTIVE:
		if (LoGZero(draw, param[4].data.d_float, param[3].data.d_int32,
				param[5].data.d_int32, param[6].data.d_int32, param[7].data.d_int32, param[8].data.d_int32) != -1)
			values[0].data.d_status = GIMP_PDB_SUCCESS;
		else
			values[0].data.d_status = GIMP_PDB_EXECUTION_ERROR;
		break;
	case GIMP_RUN_WITH_LAST_VALS:
		gimp_procedural_db_get_data("plug_in_LoG", &ldial);
		if (LoGZero(draw, ldial.sd, ldial.pa, ldial.which, ldial.pc1, ldial.pc2, param[8].data.d_int32) != -1)
			values[0].data.d_status = GIMP_PDB_SUCCESS;
		else
			values[0].data.d_status = GIMP_PDB_EXECUTION_ERROR;
		break;
	default:
		break;
	}
	gimp_drawable_detach(draw);
}

static gint LoGdialog() {
	GtkWidget *allowablePA, *logentry, *logType, *pc1entry, *pc2entry, *luminance;
	GtkWidget *dlg, *button, *frame, *main_vbox, *label, *_hbox;
	gchar *buf;
	gint i, run;
	gboolean not_done = TRUE;

	gimp_ui_init("logconv", FALSE);
	dlg = gimp_dialog_new("LoG Filter", "LoG Filter",
			NULL, 0,
			gimp_standard_help_func, "LOGCONV",
            GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
            GTK_STOCK_OK,     GTK_RESPONSE_OK,
            NULL);

	main_vbox = gtk_vbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), main_vbox, FALSE, TRUE, 0);
	gtk_widget_show(main_vbox);

	/*  Allowable PA  */
	allowablePA = gtk_combo_box_text_new();
	for (i = 0; i < 11; i++) {
		gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(allowablePA), AllowablePA[i]);
	}

	_hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(main_vbox), _hbox, FALSE, TRUE, 0);
	gtk_container_add(GTK_CONTAINER(_hbox), label = gtk_label_new("Allowable PA"));
	gtk_container_add(GTK_CONTAINER(_hbox), allowablePA);
	gtk_widget_show(label);
	gtk_widget_show(allowablePA);
	gtk_widget_show(_hbox);

	/* LoG type */
	logType = gtk_combo_box_text_new();
	for (i = 0; i < 3; i++) {
		gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(logType), LoGType[i]);
	}

	_hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(main_vbox), _hbox, FALSE, TRUE, 0);
	gtk_container_add(GTK_CONTAINER(_hbox), label = gtk_label_new("LoG Type"));
	gtk_container_add(GTK_CONTAINER(_hbox), logType);
	gtk_widget_show(label);
	gtk_widget_show(logType);
	gtk_widget_show(_hbox);

	/* Luminance algorithm */
	luminance = gtk_combo_box_text_new();
	for (i = 0; i < 3; i++) {
		gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(luminance), LuminanceType[i]);
	}

	_hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(main_vbox), _hbox, FALSE, TRUE, 0);
	gtk_container_add(GTK_CONTAINER(_hbox), label = gtk_label_new("Luminance"));
	gtk_container_add(GTK_CONTAINER(_hbox), luminance);
	gtk_widget_show(label);
	gtk_widget_show(luminance);
	gtk_widget_show(_hbox);

	/* Standard deviation */
	_hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(main_vbox), _hbox, FALSE, TRUE, 0);
	gtk_widget_show(_hbox);

	label = gtk_label_new("Standard deviation: ");
	gtk_box_pack_start(GTK_BOX(_hbox), label, FALSE, TRUE, 0);
	gtk_widget_show(label);

	buf = g_new(gchar, 10);
	sprintf(buf, "%f", ldial.sd);
	logentry = gtk_entry_new();
	gimp_help_set_help_data(logentry,
			"This is the standard deviation: the smaller it is, the finer details you'll get. The bigger it is, the coarser details you'll get. If this is smaller than 2.0, you must choose a big (3.0, 10.0) PA",
			"");
	gtk_box_pack_start(GTK_BOX(_hbox), logentry, FALSE, TRUE, 0);
	gtk_entry_set_text(GTK_ENTRY(logentry), buf);
	g_free(buf);
	gtk_widget_show(logentry);

	/* pc1 */
	_hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(main_vbox), _hbox, FALSE, TRUE, 0);
	gtk_widget_show(_hbox);

	label = gtk_label_new("PC1: ");
	gtk_box_pack_start(GTK_BOX(_hbox), label, FALSE, TRUE, 0);
	gtk_widget_show(label);

	buf = g_new(gchar, 10);
	sprintf(buf, "%d", ldial.pc1);
	pc1entry = gtk_entry_new();
	gimp_help_set_help_data(pc1entry,
			"If you choosed to apply a gradient, values below this threshold will be skipped",
			"");
	gtk_box_pack_start(GTK_BOX(_hbox), pc1entry, FALSE, TRUE, 0);
	gtk_entry_set_text(GTK_ENTRY(pc1entry), buf);
	g_free(buf);
	gtk_widget_show(pc1entry);

	/* pc2 */
	_hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(main_vbox), _hbox, FALSE, TRUE, 0);
	gtk_widget_show(_hbox);

	label = gtk_label_new("PC2: ");
	gtk_box_pack_start(GTK_BOX(_hbox), label, FALSE, TRUE, 0);
	gtk_widget_show(label);

	buf = g_new(gchar, 10);
	sprintf(buf, "%d", ldial.pc2);
	pc2entry = gtk_entry_new();
	gimp_help_set_help_data(pc2entry,
			"If you choosed to apply a gradient, values exceeding this threshold will be skipped",
			"");
	gtk_box_pack_start(GTK_BOX(_hbox), pc2entry, FALSE, TRUE, 0);
	gtk_entry_set_text(GTK_ENTRY(pc2entry), buf);
	g_free(buf);
	gtk_widget_show(pc2entry);

	gtk_widget_show(dlg);

	do {
		run = gimp_dialog_run(GIMP_DIALOG(dlg));
		if (run == GTK_RESPONSE_OK) {
			ldial.pa = gtk_combo_box_get_active(GTK_COMBO_BOX(allowablePA));
			ldial.which = gtk_combo_box_get_active(GTK_COMBO_BOX(logType));
			ldial.luminance = gtk_combo_box_get_active(GTK_COMBO_BOX(luminance));
			ldial.sd = atof(gtk_entry_get_text(GTK_ENTRY(logentry)));
			ldial.pc1 = atoi(gtk_entry_get_text(GTK_ENTRY(pc1entry)));
			ldial.pc2 = atoi(gtk_entry_get_text(GTK_ENTRY(pc2entry)));
			ldial.run = TRUE;
			gimp_set_data("plug_in_LoG", &ldial, sizeof(LoGDialog));
			not_done = FALSE;
		}
	} while (not_done);

	gtk_widget_destroy(dlg);

	return ldial.run;
}

/*
 * Computes the convolution of draw with a LoG of constant sigma and 
 * aliasing pa.
 */
static gdouble *
LoGConvolution(GimpDrawable * draw, gdouble sigma, gint pa, gint luminance) {
	gint i, j, k, addr, MSL, W, H, W1, H1, s2, s3, kd, MG, i1, j1, j2, gap;
	gdouble *temp, *temp2, *h1, *h2, *g, h1x, h1y, h2x, h2y, *temp1, Al, Ag, ks,
			sigmal, sigmag, *temp3, t;
	GimpPixelRgn src;
	guchar *buf;
	gdouble (*lFunc)(gint, gint, gint);

	W = gimp_drawable_width(draw->drawable_id);
	H = gimp_drawable_height(draw->drawable_id);
	gap = (gimp_drawable_type(draw->drawable_id) == GIMP_RGB_IMAGE) ? 3 : 1;
	buf = g_new(guchar, W * H * gap);
	gimp_pixel_rgn_init(&src, draw, 0, 0, W, H, FALSE, FALSE);
	gimp_pixel_rgn_get_rect(&src, buf, 0, 0, W, H);
	gimp_progress_init("Computing LoG...");
	/*
	 * 1. Operator design
	 */
	/* Cutoff constants */
	Ag = Table1[pa].Ag;
	Al = Table1[pa].Al;
	/* Decimation factor */
	kd = (gint) floor(KD(sigma, Ag, Al));
	/* Reconstruction constant */
	ks = KS(sigma, (gdouble ) kd, Al);
	/* Gaussian space constant */
	sigmag = SG(sigma, ks);
	/* LoG space constant */
	sigmal = SL(sigma, ks, kd);
	if ((isnan (sigmag) != 0) || (isnan (sigmal) != 0))
		return (gdouble *) NULL;
	/* LoG MSL width */
	MSL = oddsupport(sigmal);
	/* LoG filters computation */
	h1 = MakeFilter(sigmal, MSL, 1);
	h2 = MakeFilter(sigmal, MSL, 2);
	/* Gaussian MSL width */
	MG = 6 * (gint) sigmag;
	if (MG % 2 == 0)
		MG++;
	/* Gaussian filter computation */
	g = MakeFilter(sigmag, MG, 3);
	/* Luminance algorithm for RGB images */
	switch (luminance) {
	case 0:
		lFunc = LuminanceSRGB;
		break;
	case 1:
		lFunc = LuminancePerceptive1;
		break;
	case 2:
		lFunc = LuminancePerceptive2;
		break;
	default:
		return NULL;
	}
	/*
	 * 2. Convolution process
	 */
	/* Gaussian convolution and decimation */
	W1 = (gint) ceil((gdouble) W / (gdouble) kd);
	H1 = (gint) ceil((gdouble) H / (gdouble) kd);
	s3 = MG / 2;
	temp1 = g_new(gdouble, W1 * (H + s3 * 2));
	temp = g_new(gdouble, (W + s3 * 2) * H);
	i1 = W + s3 * 2;
	if (gap == 1) /* Gray images */
		for (i = 0; i < H; i++) {
			j1 = i * i1;
			for (j = 0; j < W; j++)
				temp[(j + s3) + j1] = ((gdouble) buf[j + i * W]);
		}
	else
		/* RGB images */
		for (i = 0; i < H; i++) {
			j1 = i * i1;
			for (j2 = 0, j = 0; j < W; j++, j2 += 3)
				/* Luminance */
				temp[(j + s3) + j1] = lFunc(buf[j2 + i * W * 3], buf[j2 + 1 + i * W * 3], buf[j2 + 2 + i * W * 3] * 0.114);
		}
	g_free(buf);
	gimp_progress_update(1.0 / 7.0);
	/* DC-padding */
	for (i = 0; i < H; i++) { /* Left */
		j1 = i * (W + s3 * 2);
		i1 = s3 + j1;
		for (j = 0; j < s3; j++)
			temp[j + j1] = temp[i1];
	}
	for (i = 0; i < H; i++) { /* Right */
		j1 = i * (W + s3 * 2);
		i1 = (W + s3 - 1) + j1;
		for (j = W + s3; j < (W + s3 * 2); j++)
			temp[j + j1] = temp[i1];
	}
	for (i = 0, i1 = s3; i < H; i++, i1++) /* Row convolution */
		for (j = s3, j1 = 0; j < W + s3; j += kd, j1++) {
			t = 0.0;
			for (k = 0; k < MG; k++)
				t += temp[(j - s3 + k) + i * (W + s3 * 2)] * g[k];
			temp1[j1 + i1 * W1] = t;
		}
	g_free(temp);
	gimp_progress_update(2.0 / 7.0);
	/* DC-padding */
	i1 = s3 * W1;
	for (i = 0; i < s3; i++) { /* Up */
		j1 = i * W1;
		for (j = 0; j < W1; j++)
			temp1[j + j1] = temp1[j + i1];
	}
	i1 = (H + s3 - 1) * W1;
	for (i = H + s3; i < (H + s3 * 2); i++) { /* Down */
		j1 = i * W1;
		for (j = 0; j < W1; j++)
			temp1[j + j1] = temp1[j + i1];
	}
	temp2 = g_new(gdouble, W1 * H1);
	for (i = s3, i1 = 0; i < H + s3; i += kd, i1++) /* Column convolution */
		for (j = 0; j < W1; j++) {
			t = 0.0;
			for (k = 0; k < MG; k++)
				t += temp1[j + (i - s3 + k) * W1] * g[k];
			temp2[j + i1 * W1] = t;
		}
	g_free(temp1);
	gimp_progress_update(3.0 / 7.0);
	/* LoG convolution */
	s2 = MSL / 2;
	temp = g_new(gdouble, (W1 + s2 * 2) * (H1 + s2 * 2));
	i1 = W1 + s2 * 2;
	for (i = 0; i < H1; i++)
		for (j = 0; j < W1; j++)
			temp[(j + s2) + (i + s2) * i1] = temp2[j + i * W1];
	gimp_progress_update(4.0 / 7.0);
	/* DC-padding */
	i1 = W1 + s2 * 2;
	for (i = s2; i < H1 + s2; i++) { /* Left */
		j1 = i * i1;
		for (j = 0; j < s2; j++)
			temp[j + j1] = temp[s2 + j1];
	}
	for (i = s2; i < H1 + s2; i++) { /* Right */
		j1 = i * i1;
		for (j = W1 + s2; j < (W1 + s2 * 2); j++)
			temp[j + j1] = temp[(W1 + s2 - 1) + j1];
	}
	g_free(temp2);
	temp1 = g_new(gdouble, W1 * (H1 + s2 * 2));
	temp3 = g_new(gdouble, W1 * (H1 + s2 * 2));
	for (i = s2, i1 = 0; i < H1 + s2; i++, i1++) { /* Row convolution */
		for (j = s2, j1 = 0; j < W1 + s2; j++, j1++) {
			addr = i1 * W1;
			h1x = h2x = 0.0;
			for (k = 0; k < MSL; k++) {
				h1x += temp[(j - s2 + k) + i * (W1 + s2 * 2)] * h1[k];
				h2x += temp[(j - s2 + k) + i * (W1 + s2 * 2)] * h2[k];
			}
			temp1[j1 + i * W1] = h1x;
			temp3[j1 + i * W1] = h2x;
		}
	}
	g_free(temp);
	gimp_progress_update(5.0 / 7.0);
	/* DC-padding */
	i1 = s2 * W1;
	for (i = 0; i < s2; i++) { /* Up */
		j1 = i * W1;
		for (j = 0; j < W1; j++) {
			temp1[j + j1] = temp1[j + i1];
			temp3[j + j1] = temp3[j + i1];
		}
	}
	i1 = (H1 + s2 - 1) * W1;
	for (i = H1 + s2; i < (H1 + s2 * 2); i++) { /* Down */
		j1 = i * W1;
		for (j = 0; j < W1; j++) {
			temp1[j + j1] = temp1[j + i1];
			temp3[j + j1] = temp3[j + i1];
		}
	}
	temp = g_new(gdouble, W1 * H1);
	for (i = s2, i1 = 0; i < H1 + s2; i++, i1++) { /* Column convolution */
		for (j = 0; j < W1; j++) {
			h1y = h2y = 0.0;
			for (k = 0; k < MSL; k++) {
				h1y += temp3[j + (i - s2 + k) * W1] * h1[k];
				h2y += temp1[j + (i - s2 + k) * W1] * h2[k];
			}
			temp[j + i1 * W1] = h1y + h2y;
		}
	}
	g_free(temp1);
	g_free(temp3);
	gimp_progress_update(6.0 / 7.0);
	/* Expansion */
	temp1 = g_new(gdouble, W * H);
	if (kd > 1)
		ZoomImage(temp, temp1, kd, W, H);
	else
		for (i = 0; i < H * W; i++)
			temp1[i] = temp[i];
	gimp_progress_update(7.0 / 7.0);
	g_free(temp);
	g_free(h2);
	g_free(h1);
	g_free(g);
	return temp1;
}

/*
 * Computes zero crossings for an image convoluted with a LoG of constant 
 * sigma and aliasing pa.
 * If gradient is specified, it performs a gradient thresholding as in
 * A. Basu et al. (1995).
 * gradient = 1 -> Roberts, = 2 -> Sobel
 */
static gint LoGZero(GimpDrawable * draw, gdouble sigma, gint pa, gint gradient, gint pc1, gint pc2, gint luminance) {
	gdouble *temp1, pix0, pix1, pix2, pix3, pix4, gradx, grady, *gradval, max;
	gint W, H, addr, i, j, gap, j1, i1;
	gdouble (*funcx)(gdouble *, gint, gint, gint), (*funcy)(gdouble *, gint,
			gint, gint);
	GimpPixelRgn src;
	guchar *buf;
	gboolean do_zero;

	gradval = (gdouble *) NULL;
	max = 0.0;
	funcx = NULL;
	funcy = NULL;
	W = gimp_drawable_width(draw->drawable_id);
	H = gimp_drawable_height(draw->drawable_id);
	gap = (gimp_drawable_type(draw->drawable_id) == GIMP_RGB_IMAGE) ? 3 : 1;
	if ((temp1 = LoGConvolution(draw, sigma, pa, luminance)) == (gdouble *) NULL) {
		g_message("LoG: Please choose a bigger PA.");
		return -1;
	}
	/*
	 * Zero crossing detection and gradient computation.
	 */
	if (gradient != 0) {
		max = 0.0;
		gradval = g_new(gdouble, W * H);
		memset(gradval, 0, sizeof(gdouble) * W * H);
		if (gradient == 1) {
			funcx = RobertsX;
			funcy = RobertsY;
		} else {
			funcx = SobelX;
			funcy = SobelY;
		}
	}
	buf = g_new(guchar, W * H * gap);
	memset(buf, 255, sizeof(guchar) * W * H * gap);
	gimp_progress_init("Zero crossing...");
	for (i = 1; i < H - 1; i++) {
		for (j1 = 3, j = 1; j < W - 1; j++, j1 += 3) {
			do_zero = FALSE;
			addr = j + i * W;
			/*
			 * Chessboard metric
			 *     1
			 *   2 0 3
			 *     4
			 */
			pix0 = temp1[addr];
			pix1 = temp1[j + (i - 1) * W];
			pix2 = temp1[j - 1 + i * W];
			pix3 = temp1[j + 1 + i * W];
			pix4 = temp1[j + (i + 1) * W];
			/*
			 * Zero crossing test by Simon A. J. Winder (c) 1994
			 */
			if (pix0 > 0.0
					&& (pix1 < 0.0 || pix2 < 0.0 || pix3 < 0.0 || pix4 < 0.0))
				do_zero = TRUE;
			if (pix0 == 0.0) {
				if ((pix1 > 0.0 && pix4 < 0.0) || (pix1 < 0.0 && pix4 > 0.0)
						|| (pix2 > 0.0 && pix3 < 0.0)
						|| (pix2 < 0.0 && pix3 > 0.0))
					do_zero = TRUE;
				else {
					pix1 = temp1[j - 1 + (i + 1) * W];
					pix2 = temp1[j + 1 + (i + 1) * W];
					pix3 = temp1[j - 1 + (i - 1) * W];
					pix4 = temp1[j + 1 + (i - 1) * W];
					if ((pix1 > 0.0 && pix4 < 0.0) || (pix1 < 0.0 && pix4 > 0.0)
							|| (pix2 > 0.0 && pix3 < 0.0)
							|| (pix2 < 0.0 && pix3 > 0.0))
						do_zero = TRUE;
				}
			}
			if (do_zero) {
				if (gradient != 0) {
					gradx = funcx(temp1, j, i, W);
					grady = funcy(temp1, j, i, W);
					gradval[addr] = fabs(gradx) + fabs(grady);
					if (max < gradval[addr])
						max = gradval[addr];
				} else {
					if (gap == 3)
						buf[j1 + i * W * 3] = buf[j1 + 1 + i * W * 3] = buf[j1
								+ 2 + i * W * 3] = 0;
					else
						buf[addr] = 0;
				}
			}
		}
		gimp_progress_update((gdouble) i / (gdouble) (H - 1));
	}
	if (gradient != 0) {
		/*
		 * Gradient weighting
		 */
		pc1 = (max * (gdouble) pc1) / 100.0; /* 25% */
		pc2 = (max * (gdouble) pc2) / 100.0; /* 60% */
		/*
		 * We skip all pixels below the 25% gradient threshold
		 */
		if (gap == 1) /* Gray images */
		{
			for (i = 0; i < W * H; i++)
				if ((gradval[i] > pc1) && (gradval[i] < pc2))
					buf[i] = 0;
		} else
			/* RGB images */
			for (i1 = 0, i = 0; i < W * H; i++, i1 += 3)
				if ((gradval[i] > pc1) && (gradval[i] < pc2))
					buf[i1] = buf[i1 + 1] = buf[i1 + 2] = 0;
		g_free(gradval);
	}
	gimp_pixel_rgn_init(&src, draw, 0, 0, W, H, TRUE, TRUE);
	gimp_pixel_rgn_set_rect(&src, buf, 0, 0, W, H);
	gimp_drawable_flush(draw);
	gimp_drawable_merge_shadow(draw->drawable_id, TRUE);
	gimp_drawable_update(draw->drawable_id, 0, 0, W, H);
	g_free(buf);
	g_free(temp1);
	return 0;
}

/*
 * Computes the filter's support. It must be odd.
 */
static gint oddsupport(gdouble sigma) {
	gint support;

	support = 4 * 2 * M_SQRT2 * sigma;
	if ((support % 2) == 0)
		return (support + 1);
	else
		return support;
}

/*
 * Creates filters for LoG decomposition
 */
static gdouble *
MakeFilter(gdouble sigma, gint support, gint which) {
	gint x;
	gdouble *filter;
	gdouble (*func)(gdouble, gdouble);

	func = NULL;
	switch (which) {
	case 1:
		func = H1;
		break;
	case 2:
		func = H2;
		break;
	case 3:
		func = grad;
		break;
	}
	filter = g_new(gdouble, support);
	for (x = 0; x <= support / 2; x++)
		filter[support / 2 + x] = filter[support / 2 - x] = func((gdouble) x,
				sigma);
	return filter;
}

static inline gdouble SobelX(gdouble * temp1, gint j, gint i, gint W) {
	return (-temp1[(j - 1) + (i - 1) * W] - 2.0 * temp1[j + (i - 1) * W]
			- temp1[(j + 1) + (i - 1) * W] + temp1[(j - 1) + (i + 1) * W]
			+ 2.0 * temp1[j + (i + 1) * W] + temp1[(j + 1) + (i + 1) * W]);
}

static inline gdouble SobelY(gdouble * temp1, gint j, gint i, gint W) {
	return (-temp1[(j - 1) + (i - 1) * W] - 2.0 * temp1[(j - 1) + i * W]
			- temp1[(j - 1) + (i + 1) * W] + temp1[(j + 1) + (i - 1) * W]
			+ 2.0 * temp1[(j + 1) + i * W] + temp1[(j + 1) + (i + 1) * W]);
}

static inline gdouble RobertsX(gdouble * temp1, gint j, gint i, gint W) {
	return (temp1[j + i * W] - temp1[(j + 1) + (i + 1) * W]);
}

static inline gdouble RobertsY(gdouble * temp1, gint j, gint i, gint W) {
	return (temp1[(j + 1) + i * W] - temp1[j + (i + 1) * W]);
}

static inline gdouble H1(gdouble xi, gdouble sigma) {
	return (1.0 / (sqrt(2.0 * M_PI) * pow(sigma, 2.0)))
			* (1.0 - (pow(xi, 2.0) / pow(sigma, 2.0)))
			* exp(-pow(xi, 2.0) / (2.0 * pow(sigma, 2.0)));
}

static inline gdouble H2(gdouble xi, gdouble sigma) {
	return (1.0 / (sqrt(2.0 * M_PI) * pow(sigma, 2.0)))
			* exp(-pow(xi, 2.0) / (2.0 * pow(sigma, 2.0)));
}

static inline gdouble grad(gdouble x, gdouble sigma) {
	return ((1.0 / (sqrt(2.0 * M_PI) * sigma))
			* exp(-pow(x, 2.0) / (2 * pow(sigma, 2.0))));
}

/* From http://www.w3.org/TR/AERT#color-contrast */
static inline gdouble LuminancePerceptive1(gint r, gint g, gint b) {
	return (gdouble) (r * 0.299	+ g * 0.587	+ b * 0.114);
}

/* From http://alienryderflex.com/hsp.html */
static inline gdouble LuminancePerceptive2(gint r, gint g, gint b) {
	return (gdouble) (sqrt(pow(r, 2.0) * 0.299 + pow(g, 2.0) * 0.587 + pow(b, 2.0) * 0.114));
}

/* From https://en.wikipedia.org/wiki/Luminance_%28relative%29 */
static inline gdouble LuminanceSRGB(gint r, gint g, gint b) {
    return (gdouble) (r * 0.2126 + g * 0.7152 + b * 0.0722);
}
/*
 * A C-port from LAPACK...
 */
static inline void dgefac(gdouble a[4][4], gdouble b[4]) {
	gdouble t, dmax, ipvt[4];
	gint k, i, l, j;

	/*
	 * Gaussian elimination with partial pivoting
	 */
	for (k = 0; k < 3; k++) {
		/*
		 * Look for pivot l
		 */
		l = k;
		dmax = fabs(a[k][k]);
		for (i = k + 1; i < 4; i++)
			if (fabs(a[i][k]) > dmax) {
				dmax = fabs(a[i][k]);
				l = i;
			}
		ipvt[k] = l;
		/*
		 * Swap if needed
		 */
		if (l != k) {
			t = a[l][k];
			a[l][k] = a[k][k];
			a[k][k] = t;
		}
		/*
		 * Multipliers computation
		 */
		t = -1.0 / a[k][k];
		for (i = k + 1; i < 4; i++)
			a[i][k] *= t;
		/*
		 * Row elimination
		 */
		for (j = k + 1; j < 4; j++) {
			t = a[l][j];
			if (l != k) {
				a[l][j] = a[k][j];
				a[k][j] = t;
			}
			for (i = k + 1; i < 4; i++)
				a[i][j] += (t * a[i][k]);
		}
	}
	ipvt[3] = 3;
	/*
	 * Solving L y = b
	 */
	for (k = 0; k < 3; k++) {
		l = ipvt[k];
		t = b[l];
		if (l != k) {
			b[l] = b[k];
			b[k] = t;
		}
		for (i = k + 1; i < 4; i++)
			b[i] += (t * a[i][k]);
	}
	/*
	 * Then U x = y
	 */
	for (k = 3; k >= 0; k--) {
		b[k] /= a[k][k];
		t = -b[k];
		for (i = 0; i < k; i++)
			b[i] += (a[i][k] * t);
	}
}

/*
 * Upsampling through bilinear interpolation.
 */
static void ZoomImage(gdouble * src, gdouble * dest, gint DF, gint width,
		gint height) {
	gdouble b[4], *temp;
	gdouble a[4][4];
	gint i, j, l, k, W, H, addr, addr2;

	/*
	 * Buffer initialization
	 */
	temp = g_new(gdouble, (width + DF) * (height + DF));
	memset(temp, 0, sizeof(gdouble) * (width + DF) * (height + DF));
	W = (gint) ceil((gdouble) width / (gdouble) DF);
	H = (gint) ceil((gdouble) height / (gdouble) DF);
	for (i = 0, k = 0; k < H; i += DF, k++)
		for (j = 0, l = 0; l < W; j += DF, l++)
			temp[j + i * (width + DF)] = src[l + k * W];
	/*
	 * DC-padding
	 */
	k = W * DF;
	l = (W - 1) * DF;
	for (i = 0; i < height + DF; i += DF)
		temp[k + i * (width + DF)] = temp[l + i * (width + DF)];
	k = H * DF * (width + DF);
	l = (H - 1) * DF * (width + DF);
	for (j = 0; j < width + DF; j += DF)
		temp[j + k] = temp[j + l];
	W = width + DF;
	/*
	 * Row bilinear interpolation
	 */
	for (i = 0; i < height; i += DF) {
		addr = i * W;
		addr2 = (i + DF) * W;
		for (j = 1; j < width; j += DF)
			for (k = j, l = DF - 1; k < j + DF - 1; k++, l--) {
				a[0][3] = a[1][3] = a[2][3] = a[3][3] = 1.0;
				a[0][0] = (gdouble) (k - 1);
				a[1][0] = (gdouble) (k + l);
				a[2][0] = (gdouble) (j - 1);
				a[3][0] = (gdouble) (j + DF - 1);
				a[0][1] = a[1][1] = (gdouble) i;
				a[2][1] = a[3][1] = (gdouble) (i + DF);
				a[0][2] = a[0][0] * a[0][1];
				a[1][2] = a[1][0] * a[1][1];
				a[2][2] = a[2][0] * a[2][1];
				a[3][2] = a[3][0] * a[3][1];
				b[0] = temp[(k - 1) + addr];
				b[1] = temp[(k + l) + addr];
				b[2] = temp[(j - 1) + addr2];
				b[3] = temp[(j + DF - 1) + addr2];
				dgefac(a, b);
				temp[k + addr] = (b[0] * (gdouble) k + b[1] * (gdouble) i
						+ b[2] * (gdouble) (i * k) + b[3]);
			}
	}
	/*
	 * Column bilinear interpolation
	 */
	for (i = 1; i < height; i += DF)
		for (j = 0; j < width; j++)
			for (k = i, l = DF - 1; k < i + DF - 1; k++, l--) {
				a[0][3] = a[1][3] = a[2][3] = a[3][3] = 1.0;
				a[0][0] = a[2][0] = (gdouble) j;
				a[1][0] = a[3][0] = (gdouble) (j + 1);
				a[0][1] = a[1][1] = (gdouble) (k - 1);
				a[2][1] = a[3][1] = (gdouble) (k + l);
				a[0][2] = a[0][0] * a[0][1];
				a[1][2] = a[1][0] * a[1][1];
				a[2][2] = a[2][0] * a[2][1];
				a[3][2] = a[3][0] * a[3][1];
				b[0] = temp[j + (k - 1) * W];
				b[1] = temp[(j + 1) + (k - 1) * W];
				b[2] = temp[j + (k + l) * W];
				b[3] = temp[(j + 1) + (k + l) * W];
				dgefac(a, b);
				temp[j + k * W] = (b[0] * (gdouble) j + b[1] * (gdouble) k
						+ b[2] * (gdouble) (k * j) + b[3]);
			}
	for (i = 0; i < height; i++) {
		addr = i * W;
		for (j = 0; j < width; j++)
			dest[j + i * width] = temp[j + addr];
	}
	g_free(temp);
}
