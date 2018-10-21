using System;
using System.Linq;
using System.Drawing;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using PaintDotNet;
using PaintDotNet.Effects;
using PaintDotNet.IndirectUI;
using PaintDotNet.PropertySystem;

using static System.Math;

using IntSliderControl = System.Int32;
using TextboxControl = System.String;
using ListBoxControl = System.Byte;
using DoubleSliderControl = System.Double;

[assembly: AssemblyTitle("LoGConv plugin for paint.net")]
[assembly: AssemblyDescription("Edge detection with Laplacian of Gaussian convolution")]
[assembly: AssemblyConfiguration("log")]
[assembly: AssemblyCompany("Alessandro Baldoni")]
[assembly: AssemblyProduct("LoGConv")]
[assembly: AssemblyCopyright("Copyright ©2018 by Alessandro Baldoni")]
[assembly: AssemblyTrademark("")]
[assembly: AssemblyCulture("")]
[assembly: ComVisible(false)]
[assembly: AssemblyVersion("1.0.*")]

namespace LoGConvEffect
{
    public class PluginSupportInfo : IPluginSupportInfo
    {
        public string Author
        {
            get
            {
                return ((AssemblyCopyrightAttribute)base.GetType().Assembly.GetCustomAttributes(typeof(AssemblyCopyrightAttribute), false)[0]).Copyright;
            }
        }
        public string Copyright
        {
            get
            {
                return ((AssemblyDescriptionAttribute)base.GetType().Assembly.GetCustomAttributes(typeof(AssemblyDescriptionAttribute), false)[0]).Description;
            }
        }

        public string DisplayName
        {
            get
            {
                return ((AssemblyProductAttribute)base.GetType().Assembly.GetCustomAttributes(typeof(AssemblyProductAttribute), false)[0]).Product;
            }
        }

        public Version Version
        {
            get
            {
                return base.GetType().Assembly.GetName().Version;
            }
        }

        public Uri WebsiteUri
        {
            get
            {
                return new Uri("http://www.getpaint.net/redirect/plugins.html");
            }
        }
    }

    [PluginSupportInfo(typeof(PluginSupportInfo), DisplayName = "LoG")]
    public class LoGConvEffectPlugin : PropertyBasedEffect
    {
        public static string StaticName
        {
            get
            {
                return "LoG";
            }
        }

        public static Image StaticIcon
        {
            get
            {
                return null;
            }
        }

        public static string SubmenuName
        {
            get
            {
                return "Edge Detection";
            }
        }

        public LoGConvEffectPlugin()
            : base(StaticName, StaticIcon, SubmenuName, EffectFlags.Configurable)
        {
        }

        public enum PropertyNames
        {
            AllowablePA,
            LoGType,
            StandardDeviation,
            PC1,
            PC2,
            Luminance
        }

        public enum AllowablePAOptions
        {
            Amount1Option1,
            Amount1Option2,
            Amount1Option3,
            Amount1Option4,
            Amount1Option5,
            Amount1Option6,
            Amount1Option7,
            Amount1Option8,
            Amount1Option9,
            Amount1Option10,
            Amount1Option11
        }

        public enum LoGTypeOptions
        {
            Standard,
            Roberts,
            Sobel
        }

        public enum LuminanceOptions
        {
            sRGB,
            Perceived1,
            Perceived2
        }


        protected override PropertyCollection OnCreatePropertyCollection()
        {
            List<Property> props = new List<Property>();

            AllowablePAOptions AllowablePADefault = (Enum.IsDefined(typeof(AllowablePAOptions), 0)) ? (AllowablePAOptions)0 : 0;
            props.Add(StaticListChoiceProperty.CreateForEnum<AllowablePAOptions>(PropertyNames.AllowablePA, AllowablePADefault, false));
            LoGTypeOptions LoGTypeDefault = (Enum.IsDefined(typeof(LoGTypeOptions), 0)) ? (LoGTypeOptions)0 : 0;
            props.Add(StaticListChoiceProperty.CreateForEnum<LoGTypeOptions>(PropertyNames.LoGType, LoGTypeDefault, false));
            LuminanceOptions LuminanceDefault = (Enum.IsDefined(typeof(LuminanceOptions), 0)) ? (LuminanceOptions)0 : 0;
            props.Add(StaticListChoiceProperty.CreateForEnum<LuminanceOptions>(PropertyNames.Luminance, LuminanceDefault, false));
            props.Add(new DoubleProperty(PropertyNames.StandardDeviation, 2.0, 0.0, 10.0));
            props.Add(new Int32Property(PropertyNames.PC1, 25, 0, 100));
            props.Add(new Int32Property(PropertyNames.PC2, 60, 0, 100));

            List<PropertyCollectionRule> propRules = new List<PropertyCollectionRule>();
            propRules.Add(new ReadOnlyBoundToValueRule<object, StaticListChoiceProperty>("PC1", "LoGType", LoGTypeOptions.Sobel, true));
            propRules.Add(new ReadOnlyBoundToValueRule<object, StaticListChoiceProperty>("PC2", "LoGType", LoGTypeOptions.Sobel, true));
            propRules.Add(new ReadOnlyBoundToValueRule<object, StaticListChoiceProperty>("PC1", "LoGType", LoGTypeOptions.Roberts, true));
            propRules.Add(new ReadOnlyBoundToValueRule<object, StaticListChoiceProperty>("PC2", "LoGType", LoGTypeOptions.Roberts, true));

            return new PropertyCollection(props, propRules);
        }

        protected override ControlInfo OnCreateConfigUI(PropertyCollection props)
        {
            ControlInfo configUI = CreateDefaultConfigUI(props);

            configUI.SetPropertyControlValue(PropertyNames.AllowablePA, ControlInfoPropertyNames.DisplayName, "Allowable PA");
            PropertyControlInfo Amount1Control = configUI.FindControlForPropertyName(PropertyNames.AllowablePA);
            Amount1Control.SetValueDisplayName(AllowablePAOptions.Amount1Option1, "0.0001");
            Amount1Control.SetValueDisplayName(AllowablePAOptions.Amount1Option2, "0.0003");
            Amount1Control.SetValueDisplayName(AllowablePAOptions.Amount1Option3, "0.0010");
            Amount1Control.SetValueDisplayName(AllowablePAOptions.Amount1Option4, "0.0030");
            Amount1Control.SetValueDisplayName(AllowablePAOptions.Amount1Option5, "0.0100");
            Amount1Control.SetValueDisplayName(AllowablePAOptions.Amount1Option6, "0.0250");
            Amount1Control.SetValueDisplayName(AllowablePAOptions.Amount1Option7, "0.1000");
            Amount1Control.SetValueDisplayName(AllowablePAOptions.Amount1Option8, "0.3000");
            Amount1Control.SetValueDisplayName(AllowablePAOptions.Amount1Option9, "1.0000");
            Amount1Control.SetValueDisplayName(AllowablePAOptions.Amount1Option10, "3.0000");
            Amount1Control.SetValueDisplayName(AllowablePAOptions.Amount1Option11, "10.0000");
            configUI.SetPropertyControlValue(PropertyNames.LoGType, ControlInfoPropertyNames.DisplayName, "LoG Type");
            PropertyControlInfo Amount2Control = configUI.FindControlForPropertyName(PropertyNames.LoGType);
            Amount2Control.SetValueDisplayName(LoGTypeOptions.Standard, "Standard LoG");
            Amount2Control.SetValueDisplayName(LoGTypeOptions.Roberts, "LoG with Roberts");
            Amount2Control.SetValueDisplayName(LoGTypeOptions.Sobel, "LoG with Sobel");
            PropertyControlInfo LuminanceControl = configUI.FindControlForPropertyName(PropertyNames.Luminance);
            LuminanceControl.SetValueDisplayName(LuminanceOptions.sRGB, "sRGB");
            LuminanceControl.SetValueDisplayName(LuminanceOptions.Perceived1, "Perceived 1");
            LuminanceControl.SetValueDisplayName(LuminanceOptions.Perceived2, "Perceived 2");
            configUI.SetPropertyControlValue(PropertyNames.StandardDeviation, ControlInfoPropertyNames.DisplayName, "Standard Deviation");
            configUI.SetPropertyControlValue(PropertyNames.PC1, ControlInfoPropertyNames.DisplayName, "PC1");
            configUI.SetPropertyControlValue(PropertyNames.PC2, ControlInfoPropertyNames.DisplayName, "PC2");

            return configUI;
        }

        protected override void OnCustomizeConfigUIWindowProperties(PropertyCollection props)
        {
            // Change the effect's window title
            props[ControlInfoPropertyNames.WindowTitle].Value = "LoG Parameters";
            // Add help button to effect UI
            props[ControlInfoPropertyNames.WindowHelpContentType].Value = WindowHelpContentType.PlainText;
            props[ControlInfoPropertyNames.WindowHelpContent].Value = " v1,0\nCopyright ©2018 by \nAll rights reserved.";
            base.OnCustomizeConfigUIWindowProperties(props);
        }

        protected override void OnSetRenderInfo(PropertyBasedEffectConfigToken newToken, RenderArgs dstArgs, RenderArgs srcArgs)
        {
            AllowablePA = (byte)((int)newToken.GetProperty<StaticListChoiceProperty>(PropertyNames.AllowablePA).Value);
            LoGType = (byte)((int)newToken.GetProperty<StaticListChoiceProperty>(PropertyNames.LoGType).Value);
            StandardDeviation = newToken.GetProperty<DoubleProperty>(PropertyNames.StandardDeviation).Value;
            PC1 = newToken.GetProperty<Int32Property>(PropertyNames.PC1).Value;
            PC2 = newToken.GetProperty<Int32Property>(PropertyNames.PC2).Value;
            Luminance = (byte)((int)newToken.GetProperty<StaticListChoiceProperty>(PropertyNames.Luminance).Value);

            base.OnSetRenderInfo(newToken, dstArgs, srcArgs);
        }

        protected override unsafe void OnRender(Rectangle[] rois, int startIndex, int length)
        {
            if (length == 0) return;
            for (int i = startIndex; i < startIndex + length; ++i)
            {
                Render(DstArgs.Surface, SrcArgs.Surface, rois[i]);
            }
        }

        #region User Entered Code
        // Name:
        // Submenu:
        // Author:
        // Title:
        // Version:
        // Desc:
        // Keywords:
        // URL:
        // Help:
        #region UICode
        ListBoxControl AllowablePA = 0; // Allowable PA|0.0001|0.0003|0.0010|0.0030|0.0100|0.0250|0.1000|0.3000|1.0000|3.0000|10.0000
        ListBoxControl LoGType = 0; // LoG Type|Standard LoG|LoG with Roberts|LoG with Sobel
        DoubleSliderControl StandardDeviation = 2.0; // [0,255] Standard Deviation
        IntSliderControl PC1 = 25; // [0,100] PC1
        IntSliderControl PC2 = 60; // [0,100] PC2
        ListBoxControl Luminance = 0; //Luminance|sRGB|Perceived 1|Perceived 2
        #endregion

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double KD(double sigma, double ag, double al)
        {
            return ((sigma * PI) / (Sqrt(Pow(ag, 2.0) + Pow(al, 2.0))));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double KS(double sigma, double kd, double al)
        {
            return ((sigma * PI) / (al * kd));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double SG(double sigma, double ks)
        {
            return (sigma * Sqrt(1.0 - (1.0 / Pow(ks, 2.0))));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double SL(double sigma, double ks, double kd)
        {
            return (sigma / (ks * kd));
        }

        private static readonly ValueTuple<double, double>[] Table1 = new ValueTuple<double, double>[]{
            new ValueTuple<double, double>(3.554140, 4.205983),		// pa =  0.0001 
            new ValueTuple<double, double>(3.402057, 4.061595),		// pa =  0.0003 
            new ValueTuple<double, double>(3.227792, 3.896178),		// pa =  0.0010 
            new ValueTuple<double, double>(3.060844, 3.737694),		// pa =  0.0030 
            new ValueTuple<double, double>(2.867757, 3.554300),		// pa =  0.0100 
            new ValueTuple<double, double>(2.712526, 3.406712),		// pa =  0.0250 
            new ValueTuple<double, double>(2.461219, 3.167270),		// pa =  0.1000 
            new ValueTuple<double, double>(2.244686, 2.960158),		// pa =  0.3000 
            new ValueTuple<double, double>(1.984301, 2.709568),		// pa =  1.0000 
            new ValueTuple<double, double>(1.718008, 2.450782),		// pa =  3.0000 
            new ValueTuple<double, double>(1.378026, 2.115151)		// pa = 10.0000 
         };

        private static readonly double M_SQRT2 = Sqrt(2.0);

        /*
         * Computes the filter's support. It must be odd.
         */
        private static int oddsupport(double sigma)
        {
            int support;

            support = (int)(4.0 * 2.0 * M_SQRT2 * sigma);
            if ((support % 2) == 0)
                return (support + 1);
            else
                return support;
        }

        /*
         * Creates filters for LoG decomposition
         */
        private static double[] MakeFilter(double sigma, int support, int which)
        {
            int x;
            double[] filter;
            Func<double, double, double> func;

            func = null;
            switch (which)
            {
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
            filter = new double[support];
            for (x = 0; x <= support / 2; x++)
                filter[support / 2 + x] = filter[support / 2 - x] = func((double)x, sigma);
            return filter;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double SobelX(double[] temp1, int j, int i, int W)
        {
            return (-temp1[(j - 1) + (i - 1) * W] - 2.0 * temp1[j + (i - 1) * W] -
                temp1[(j + 1) + (i - 1) * W] + temp1[(j - 1) + (i + 1) * W] +
                2.0 * temp1[j + (i + 1) * W] + temp1[(j + 1) + (i + 1) * W]);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double SobelY(double[] temp1, int j, int i, int W)
        {
            return (-temp1[(j - 1) + (i - 1) * W] - 2.0 * temp1[(j - 1) + i * W] -
                temp1[(j - 1) + (i + 1) * W] + temp1[(j + 1) + (i - 1) * W] +
                2.0 * temp1[(j + 1) + i * W] + temp1[(j + 1) + (i + 1) * W]);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double RobertsX(double[] temp1, int j, int i, int W)
        {
            return (temp1[j + i * W] - temp1[(j + 1) + (i + 1) * W]);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double RobertsY(double[] temp1, int j, int i, int W)
        {
            return (temp1[(j + 1) + i * W] - temp1[j + (i + 1) * W]);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double H1(double xi, double sigma)
        {
            return (1.0 / (Sqrt(2.0 * PI) * Pow(sigma, 2.0))) *
              (1.0 - (Pow(xi, 2.0) / Pow(sigma, 2.0))) *
              Exp(-Pow(xi, 2.0) / (2.0 * Pow(sigma, 2.0)));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double H2(double xi, double sigma)
        {
            return (1.0 / (Sqrt(2.0 * PI) * Pow(sigma, 2.0))) *
              Exp(-Pow(xi, 2.0) / (2.0 * Pow(sigma, 2.0)));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double grad(double x, double sigma)
        {
            return ((1.0 / (Sqrt(2.0 * PI) * sigma)) *
                Exp(-Pow(x, 2.0) / (2 * Pow(sigma, 2.0))));
        }

        // From https://en.wikipedia.org/wiki/Luminance_%28relative%29
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double luminance_sRGB(byte r, byte g, byte b)
        {
            return ((double)(r * 0.2126 + g * 0.7152 + b * 0.0722));
        }

        // From http://www.w3.org/TR/AERT#color-contrast
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double luminance_Perceived1(byte r, byte g, byte b)
        {
            return ((double)(r * 0.299 + g * 0.587 + b * 0.114));
        }

        // From http://alienryderflex.com/hsp.html
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double luminance_Perceived2(byte r, byte g, byte b)
        {
            return (Sqrt(Pow(r, 2.0) * 0.299 + Pow(g, 2.0) * 0.587 + Pow(b, 2.0) * 0.114));
        }

        /*
         * A C-port from LAPACK...
         */
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void dgefac(double[,] a, double[] b)
        {
            double t, dmax;
            int[] ipvt = new int[4];
            int k, i, l, j;

            /*
             * Gaussian elimination with partial pivoting 
             */
            for (k = 0; k < 3; k++)
            {
                /*
                 * Look for pivot l
                 */
                l = k;
                dmax = Abs(a[k, k]);
                for (i = k + 1; i < 4; i++)
                    if (Abs(a[i, k]) > dmax)
                    {
                        dmax = Abs(a[i, k]);
                        l = i;
                    }
                ipvt[k] = l;
                /*
                 * Swap if needed
                 */
                if (l != k)
                {
                    t = a[l, k];
                    a[l, k] = a[k, k];
                    a[k, k] = t;
                }
                /*
                 * Multipliers computation
                 */
                t = -1.0 / a[k, k];
                for (i = k + 1; i < 4; i++)
                    a[i, k] *= t;
                /*
                 * Row elimination
                 */
                for (j = k + 1; j < 4; j++)
                {
                    t = a[l, j];
                    if (l != k)
                    {
                        a[l, j] = a[k, j];
                        a[k, j] = t;
                    }
                    for (i = k + 1; i < 4; i++)
                        a[i, j] += (t * a[i, k]);
                }
            }
            ipvt[3] = 3;
            /*
             * Solving L y = b
             */
            for (k = 0; k < 3; k++)
            {
                l = ipvt[k];
                t = b[l];
                if (l != k)
                {
                    b[l] = b[k];
                    b[k] = t;
                }
                for (i = k + 1; i < 4; i++)
                    b[i] += (t * a[i, k]);
            }
            /*
             * Then U x = y
             */
            for (k = 3; k >= 0; k--)
            {
                b[k] /= a[k, k];
                t = -b[k];
                for (i = 0; i < k; i++)
                    b[i] += (a[i, k] * t);
            }
        }


        /*
         * Upsampling through bilinear interpolation.
         */
        static void
        ZoomImage(double[] src, double[] dest, int DF, int width, int height)
        {
            double[] b = new double[4];
            double[] temp;
            double[,] a = new double[4, 4];
            int i, j, l, k, W, H, addr, addr2;

            /*
             * Buffer initialization
             */
            temp = new double[(width + DF) * (height + DF)];
            //memset(temp, 0, sizeof(double) * (width + DF) * (height + DF));
            W = (int)Ceiling((double)width / (double)DF);
            H = (int)Ceiling((double)height / (double)DF);
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
            for (i = 0; i < height; i += DF)
            {
                addr = i * W;
                addr2 = (i + DF) * W;
                for (j = 1; j < width; j += DF)
                    for (k = j, l = DF - 1; k < j + DF - 1; k++, l--)
                    {
                        a[0, 3] = a[1, 3] = a[2, 3] = a[3, 3] = 1.0;
                        a[0, 0] = (double)(k - 1);
                        a[1, 0] = (double)(k + l);
                        a[2, 0] = (double)(j - 1);
                        a[3, 0] = (double)(j + DF - 1);
                        a[0, 1] = a[1, 1] = (double)i;
                        a[2, 1] = a[3, 1] = (double)(i + DF);
                        a[0, 2] = a[0, 0] * a[0, 1];
                        a[1, 2] = a[1, 0] * a[1, 1];
                        a[2, 2] = a[2, 0] * a[2, 1];
                        a[3, 2] = a[3, 0] * a[3, 1];
                        b[0] = temp[(k - 1) + addr];
                        b[1] = temp[(k + l) + addr];
                        b[2] = temp[(j - 1) + addr2];
                        b[3] = temp[(j + DF - 1) + addr2];
                        dgefac(a, b);
                        temp[k + addr] = (b[0] * (double)k + b[1] * (double)i + b[2] * (double)(i * k) + b[3]);
                    }
            }
            /*
             * Column bilinear interpolation 
             */
            for (i = 1; i < height; i += DF)
                for (j = 0; j < width; j++)
                    for (k = i, l = DF - 1; k < i + DF - 1; k++, l--)
                    {
                        a[0, 3] = a[1, 3] = a[2, 3] = a[3, 3] = 1.0;
                        a[0, 0] = a[2, 0] = (double)j;
                        a[1, 0] = a[3, 0] = (double)(j + 1);
                        a[0, 1] = a[1, 1] = (double)(k - 1);
                        a[2, 1] = a[3, 1] = (double)(k + l);
                        a[0, 2] = a[0, 0] * a[0, 1];
                        a[1, 2] = a[1, 0] * a[1, 1];
                        a[2, 2] = a[2, 0] * a[2, 1];
                        a[3, 2] = a[3, 0] * a[3, 1];
                        b[0] = temp[j + (k - 1) * W];
                        b[1] = temp[(j + 1) + (k - 1) * W];
                        b[2] = temp[j + (k + l) * W];
                        b[3] = temp[(j + 1) + (k + l) * W];
                        dgefac(a, b);
                        temp[j + k * W] = (b[0] * (double)j + b[1] * (double)k +
                                         b[2] * (double)(k * j) + b[3]);
                    }
            for (i = 0; i < height; i++)
            {
                addr = i * W;
                for (j = 0; j < width; j++)
                    dest[j + i * width] = temp[j + addr];
            }
        }

        /*
         * Computes the convolution of draw with a LoG of constant sigma and 
         * aliasing pa.
         */
        private double[] LoGConvolution(Surface src, double sigma, int pa, int luminanceAlgorithm)
        {
            int i, j, k, addr, MSL, W, H, W1, H1, s2, s3, kd, MG, i1, j1, j2, gap;
            double[] temp, temp2, h1, h2, g, temp1, temp3;
            double h1x, h1y, h2x, h2y, Al, Ag, ks, sigmal, sigmag, t;
            byte[] buf;
            Func<byte, byte, byte, double> luminance = null;

            W = src.Width;
            H = src.Height;
            // Original code:
            //   gap = (gimp_drawable_type(draw->id) == RGB_IMAGE) ? 3 : 1;
            // In Paint.NET, all images are RGBa
            gap = 4;
            buf = src.Scan0.ToByteArray();
            /* 
             * 1. Operator design 
             */
            /* Cutoff constants */
            Ag = Table1[pa].Item1;
            Al = Table1[pa].Item2;
            /* Decimation factor */
            kd = (int)Floor(KD(sigma, Ag, Al));
            /* Reconstruction constant */
            ks = KS(sigma, (double)kd, Al);
            /* Gaussian space constant */
            sigmag = SG(sigma, ks);
            /* LoG space constant */
            sigmal = SL(sigma, ks, kd);
            if (Double.IsNaN(sigmag) || Double.IsNaN(sigmal))
                return null;
            /* LoG MSL width */
            MSL = oddsupport(sigmal);
            /* LoG filters computation */
            h1 = MakeFilter(sigmal, MSL, 1);
            h2 = MakeFilter(sigmal, MSL, 2);
            /* Gaussian MSL width */
            MG = 6 * (int)sigmag;
            if (MG % 2 == 0)
                MG++;
            /* Gaussian filter computation */
            g = MakeFilter(sigmag, MG, 3);
            /* Luminance algorithm */
            switch (luminanceAlgorithm)
            {
                case 0:
                    luminance = luminance_sRGB;
                    break;
                case 1:
                    luminance = luminance_Perceived1;
                    break;
                case 2:
                    luminance = luminance_Perceived2;
                    break;
            }
            /*
             * 2. Convolution process
             */
            /* Gaussian convolution and decimation */
            W1 = (int)Ceiling((double)W / (double)kd);
            H1 = (int)Ceiling((double)H / (double)kd);
            s3 = MG / 2;
            temp1 = new double[W1 * (H + s3 * 2)];
            temp = new double[(W + s3 * 2) * H];
            i1 = W + s3 * 2;
            if (gap == 1)                 /* Gray images */
                for (i = 0; i < H; i++)
                {
                    j1 = i * i1;
                    for (j = 0; j < W; j++)
                        temp[(j + s3) + j1] = ((double)buf[j + i * W]);
                }
            else                          /* RGB images */
                for (i = 0; i < H; i++)
                {
                    j1 = i * i1;
                    for (j2 = 0, j = 0; j < W; j++, j2 += 4)
                        /* Luminance */
                        // In paint.net pixels are in BGRa order
                        temp[(j + s3) + j1] = luminance(buf[j2 + 2 + i * W * 4], buf[j2 + 1 + i * W * 4], buf[j2 + 2 + i * W * 4]);
                }
            /* DC-padding */
            for (i = 0; i < H; i++)
            {               /* Left */
                j1 = i * (W + s3 * 2);
                i1 = s3 + j1;
                for (j = 0; j < s3; j++)
                    temp[j + j1] = temp[i1];
            }
            for (i = 0; i < H; i++)
            {               /* Right */
                j1 = i * (W + s3 * 2);
                i1 = (W + s3 - 1) + j1;
                for (j = W + s3; j < (W + s3 * 2); j++)
                    temp[j + j1] = temp[i1];
            }
            for (i = 0, i1 = s3; i < H; i++, i1++)  /* Row convolution */
                for (j = s3, j1 = 0; j < W + s3; j += kd, j1++)
                {
                    t = 0.0;
                    for (k = 0; k < MG; k++)
                        t += temp[(j - s3 + k) + i * (W + s3 * 2)] * g[k];
                    temp1[j1 + i1 * W1] = t;
                }
            /* DC-padding */
            i1 = s3 * W1;
            for (i = 0; i < s3; i++)
            {               /* Up */
                j1 = i * W1;
                for (j = 0; j < W1; j++)
                    temp1[j + j1] = temp1[j + i1];
            }
            i1 = (H + s3 - 1) * W1;
            for (i = H + s3; i < (H + s3 * 2); i++)
            {               /* Down */
                j1 = i * W1;
                for (j = 0; j < W1; j++)
                    temp1[j + j1] = temp1[j + i1];
            }
            temp2 = new double[W1 * H1];
            for (i = s3, i1 = 0; i < H + s3; i += kd, i1++)     /* Column convolution */
                for (j = 0; j < W1; j++)
                {
                    t = 0.0;
                    for (k = 0; k < MG; k++)
                        t += temp1[j + (i - s3 + k) * W1] * g[k];
                    temp2[j + i1 * W1] = t;
                }
            /* LoG convolution */
            s2 = MSL / 2;
            temp = new double[(W1 + s2 * 2) * (H1 + s2 * 2)];
            i1 = W1 + s2 * 2;
            for (i = 0; i < H1; i++)
                for (j = 0; j < W1; j++)
                    temp[(j + s2) + (i + s2) * i1] = temp2[j + i * W1];
            /* DC-padding */
            i1 = W1 + s2 * 2;
            for (i = s2; i < H1 + s2; i++)
            {               /* Left */
                j1 = i * i1;
                for (j = 0; j < s2; j++)
                    temp[j + j1] = temp[s2 + j1];
            }
            for (i = s2; i < H1 + s2; i++)
            {               /* Right */
                j1 = i * i1;
                for (j = W1 + s2; j < (W1 + s2 * 2); j++)
                    temp[j + j1] = temp[(W1 + s2 - 1) + j1];
            }
            temp1 = new double[W1 * (H1 + s2 * 2)];
            temp3 = new double[W1 * (H1 + s2 * 2)];
            for (i = s2, i1 = 0; i < H1 + s2; i++, i1++)
            {               /* Row convolution */
                for (j = s2, j1 = 0; j < W1 + s2; j++, j1++)
                {
                    addr = i1 * W1;
                    h1x = h2x = 0.0;
                    for (k = 0; k < MSL; k++)
                    {
                        h1x += temp[(j - s2 + k) + i * (W1 + s2 * 2)] * h1[k];
                        h2x += temp[(j - s2 + k) + i * (W1 + s2 * 2)] * h2[k];
                    }
                    temp1[j1 + i * W1] = h1x;
                    temp3[j1 + i * W1] = h2x;
                }
            }
            /* DC-padding */
            i1 = s2 * W1;
            for (i = 0; i < s2; i++)
            {               /* Up */
                j1 = i * W1;
                for (j = 0; j < W1; j++)
                {
                    temp1[j + j1] = temp1[j + i1];
                    temp3[j + j1] = temp3[j + i1];
                }
            }
            i1 = (H1 + s2 - 1) * W1;
            for (i = H1 + s2; i < (H1 + s2 * 2); i++)
            {               /* Down */
                j1 = i * W1;
                for (j = 0; j < W1; j++)
                {
                    temp1[j + j1] = temp1[j + i1];
                    temp3[j + j1] = temp3[j + i1];
                }
            }
            temp = new double[W1 * H1];
            for (i = s2, i1 = 0; i < H1 + s2; i++, i1++)
            {               /* Column convolution */
                for (j = 0; j < W1; j++)
                {
                    h1y = h2y = 0.0;
                    for (k = 0; k < MSL; k++)
                    {
                        h1y += temp3[j + (i - s2 + k) * W1] * h1[k];
                        h2y += temp1[j + (i - s2 + k) * W1] * h2[k];
                    }
                    temp[j + i1 * W1] = h1y + h2y;
                }
            }
            /* Expansion */
            temp1 = new double[W * H];
            if (kd > 1)
                ZoomImage(temp, temp1, kd, W, H);
            else
                for (i = 0; i < H * W; i++)
                    temp1[i] = temp[i];
            return temp1;
        }

        /*
         * Computes zero crossings for an image convoluted with a LoG of constant 
         * sigma and aliasing pa.
         * If gradient is specified, it performs a gradient thresholding as in
         * A. Basu et al. (1995).
         * gradient = 1 -> Roberts, = 2 -> Sobel
         */
        private int LoGZero(Surface src, Surface dest, double sigma, int pa, int gradient, int luminance)
        {
            double[] temp1, gradval;
            double pix0, pix1, pix2, pix3, pix4, gradx, grady, max, pc1, pc2;
            int W, H, addr, i, j, gap, j1, i1;
            Func<double[], int, int, int, double> funcx;
            Func<double[], int, int, int, double> funcy;
            byte[] buf;
            bool do_zero;

            gradval = null;
            max = 0.0;
            funcx = null;
            funcy = null;
            W = src.Width;
            H = src.Height;
            // Original code:
            //      gap = (gimp_drawable_type(draw->id) == RGB_IMAGE) ? 3 : 1;
            // In Paint.NET all images are RGBa
            gap = 4;
            if ((temp1 = LoGConvolution(src, sigma, pa, luminance)) == null)
            {
                //g_message("LoG: Please choose a bigger PA.");
                return -1;
            }
            /*
             * Zero crossing detection and gradient computation.
             */
            if (gradient != 0)
            {
                max = 0.0;
                gradval = new double[W * H];
                if (gradient == 1)
                {
                    funcx = RobertsX;
                    funcy = RobertsY;
                }
                else
                {
                    funcx = SobelX;
                    funcy = SobelY;
                }
            }
            buf = Enumerable.Repeat<byte>(255, W * H * gap).ToArray();
            for (i = 1; i < H - 1; i++)
            {
                for (j1 = 4, j = 1; j < W - 1; j++, j1 += 4)
                {
                    do_zero = false;
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
                    if (pix0 > 0.0 &&
                        (pix1 < 0.0 || pix2 < 0.0 || pix3 < 0.0 || pix4 < 0.0))
                        do_zero = true;
                    if (pix0 == 0.0)
                    {
                        if ((pix1 > 0.0 && pix4 < 0.0) || (pix1 < 0.0 && pix4 > 0.0) ||
                        (pix2 > 0.0 && pix3 < 0.0) || (pix2 < 0.0 && pix3 > 0.0))
                            do_zero = true;
                        else
                        {
                            pix1 = temp1[j - 1 + (i + 1) * W];
                            pix2 = temp1[j + 1 + (i + 1) * W];
                            pix3 = temp1[j - 1 + (i - 1) * W];
                            pix4 = temp1[j + 1 + (i - 1) * W];
                            if ((pix1 > 0.0 && pix4 < 0.0) ||
                                (pix1 < 0.0 && pix4 > 0.0) ||
                                (pix2 > 0.0 && pix3 < 0.0) || (pix2 < 0.0 && pix3 > 0.0))
                                do_zero = true;
                        }
                    }
                    if (do_zero)
                    {
                        if (gradient != 0)
                        {
                            gradx = funcx(temp1, j, i, W);
                            grady = funcy(temp1, j, i, W);
                            gradval[addr] = Abs(gradx) + Abs(grady);
                            if (max < gradval[addr])
                                max = gradval[addr];
                        }
                        else
                        {
                            if (gap > 1)
                            {
                                buf[j1 + i * W * 4] = buf[j1 + 1 + i * W * 4] = buf[j1 + 2 + i * W * 4] = 0;
                                // Alpha channel set to opaque
                                buf[j1 + 3 + i * W * 4] = 255;
                            }
                            else
                                buf[addr] = 0;
                        }
                    }
                }
            }
            if (gradient != 0)
            {
                /*
                 * Gradient weighting
                 */
                pc1 = (max * (double)PC1) / 100.0;   /* 25% */
                pc2 = (max * (double)PC2) / 100.0;   /* 60% */
                                                     /*
                                                      * We skip all pixels below the 25% gradient threshold
                                                      */
                if (gap == 1)          /* Gray images */
                {
                    for (i = 0; i < W * H; i++)
                        if ((gradval[i] > pc1) && (gradval[i] < pc2))
                            buf[i] = 0;
                }
                else                   /* RGBa images */
                    for (i1 = 0, i = 0; i < W * H; i++, i1 += 4)
                        if ((gradval[i] > pc1) && (gradval[i] < pc2))
                        {
                            buf[i1] = buf[i1 + 1] = buf[i1 + 2] = 0;
                            // Alpha channel is set to opaque
                            buf[i1 + 3] = 255;
                        }
            }


            unsafe
            {
                fixed (byte* bufPtr = &buf[0])
                {
                    Buffer.MemoryCopy(bufPtr, dest.Scan0.VoidStar, dest.Scan0.Length, dest.Scan0.Length);
                }
            }
            return 0;
        }


        void Render(Surface dst, Surface src, Rectangle rect)
        {
            LoGZero(src, dst, StandardDeviation, AllowablePA, LoGType, Luminance);
        }

        #endregion
    }
}
