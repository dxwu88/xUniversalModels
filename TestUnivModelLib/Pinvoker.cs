using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Microsoft.Win32.SafeHandles;
using System.Security.Permissions;
using System.Runtime.InteropServices;
using System.Runtime.ConstrainedExecution;
using System.IO;

namespace Honeywell.UOP.CPP2CS
{
    public class PInvoker
    {
        static PInvoker()
        {
            string pathToSet = Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location);
            SetDllDirectory(pathToSet);
        }

        [DllImport(@"xUniversalModels.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void RunUnivModel(TemperatureProfileEnum tpe, int ysize, double[] yin, double[] yout, int tsize, double[] temperatureProfile);

        [DllImport(@"xUniversalModels.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void TestMPI(int size, int[] a, int[] b);

        [DllImport(@"C:\working_copy\TestCPP\Debug\TestCPP.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void test1(int[] a, int[] b);

        [DllImport(@"C:\working_copy\TestCPP\Debug\TestCPP.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern bool test2(int[,] a);

        [DllImport(@"xUniversalModels.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void test3(int size, double[] a, double[] b);

        [DllImport(@"C:\working_copy\TestCPP\Debug\TestCPP.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void test4(ref UserDataStruct abc);

        [DllImport(@"C:\working_copy\TestCPP\Debug\TestCPP.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void test5(string strx, StringBuilder sb, int c);

        [DllImport(@"C:\working_copy\TestCPP\Debug\TestCPP.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern IntPtr test6(string strx);


        [DllImport("kernel32.dll", CharSet = CharSet.Unicode, SetLastError = true)]
        [return: MarshalAs(UnmanagedType.Bool)]
        public static extern bool SetDllDirectory(string lpPathName);
    }

}

