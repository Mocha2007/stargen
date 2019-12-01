using System;
using System.IO;
// for windows
using System.ComponentModel;
using System.Drawing;
using System.Windows.Forms;
// for pictures
// using System.Drawing.Image;

public class Program
{
	static Random rnd = new Random();
	public static StarSystem system = RandSystem();
	static void Main(string[] args)
	{
		Console.CursorVisible = false;
		Application.EnableVisualStyles();
		Application.Run(new Interface());
		// Console.ReadKey();
	}
	// rng
	static double RandAngle(){
		return Uniform(0, 2*Math.PI);
	}
	static double RandUnit(){
		// MAX_INT = 0x7FFFFFFF;
		return (double)rnd.Next(0, 0x7FFFFFFF) / 0x7FFFFFFF;
	}
	static double Uniform(double minimum, double maximum){
		return RandUnit() * (maximum - minimum) + minimum;
	}
	static double LogUniform(double minimum, double maximum){
		return Math.Exp(Uniform(Math.Log(minimum), Math.Log(maximum)));
	}
	// random astro - todo
	static StarSystem RandSystem(){
		Star primary = RandStar();
		byte planet_n = (byte)rnd.Next(7, 10);
		Planet[] secondaries = new Planet[planet_n];
		double sma = Math.Sqrt(primary.luminosity / Constants.sun.luminosity) * Constants.mercury.orbit.sma;
		for (byte i = 0; i < planet_n; i++){
			secondaries[i] = RandPlanet(primary, sma);
			sma *= Uniform(1.38, 2.01); // bode's law... kinda...
		}
		return new StarSystem(primary, secondaries);
	}
	static Star RandStar(){
		double mass = LogUniform(0.08, 5);
		double radius = Constants.sun.radius * Math.Pow(mass, 0.96);
		double luminosity = 0.45 < mass ? Math.Pow(1.148 * mass, 3.4751) : Math.Pow(0.2264 * mass, 2.52);
		ushort temperature = (ushort)(5772 * Math.Pow(mass, 0.54));
		return new Star(Constants.sun.mass * mass, radius, Constants.sun.luminosity * luminosity, temperature);
	}
	static Planet RandPlanet(Star primary, double sma){
		double mass, albedo, density;
		double aop = RandAngle();
		double ecc = Uniform(0, 0.21);
		double man = RandAngle();
		Orbit orbit = new Orbit(primary, aop, ecc, man, sma);
		if (primary.FrostLine() < sma){ //giant
			mass = LogUniform(10*Constants.earth.mass, 10*Constants.jupiter.mass);
			if (2e26 < mass){ // gas giant
				density = Uniform(600, 1400);
				albedo = Uniform(0.34, 0.5);
			}
			else{ // ice giant
				density = Uniform(1200, 1700);
				albedo = Uniform(0.29, 0.3);
			}
		}
		else{ // rocky
			mass = LogUniform(Constants.mercury.mass / 20, 10*Constants.earth.mass);
			density = Uniform(3900, 5600);
			albedo = Uniform(0.09, 0.76);
		}
		double radius = Math.Pow(mass/density * 3/4 / Math.PI, (double)1/3);
		return new Planet(mass, radius, albedo, orbit);
	}
}
public class Constants
{
	public static float G = 6.6743e-11F;
	public static ulong au = 149597870700;
	public static Star sun = new Star(1.9885e30, 6.957e8, 3.828e26, 5778);
	static Orbit mercuryOrbit = new Orbit(sun, 0.50831, 0.20563, 3.05077, 5.790905e10);
	public static Planet mercury = new Planet(3.3011e23, 2.4397e6, 0.088, mercuryOrbit);
	static Orbit earthOrbit = new Orbit(sun, 0, 0.0167086, 0, 1.49598023e11);
	public static Planet earth = new Planet(5.97237e24, 6.371e6, 0.306, earthOrbit);
	static Orbit jupiterOrbit = new Orbit(sun, 4.77988, 0.0489, 0.34941, 5.2044 * au);
	public static Planet jupiter = new Planet(1.8982e27, 6.9911e7, 0.503, jupiterOrbit);
	public static double Remap(double value, double[] from, double[] to){
		double from_min = from[0];
		double from_max = from[1];
		double from_range = from_max - from_min;
		double to_min = to[0];
		double to_max = to[1];
		double to_range = to_max - to_min;
		return (value - from_min)/from_range * to_range + to_min;
	}
}
public class Body
{
	public double mass, radius;
	public Body(double mass, double radius){
		this.mass = mass;
		this.radius = radius;
	}
	public override string ToString(){
		return "Body {\n\tmass = " + this.mass.ToString() + " kg\n\tradius = " + (this.radius/1000).ToString() + " km\n}";
	}
	public double Mu(){
		return Constants.G * this.mass;
	}
}
public class Planet : Body
{
	public double albedo;
	public Orbit orbit;
	public Planet(double mass, double radius, double albedo, Orbit orbit) : base(mass, radius){
		this.mass = mass;
		this.radius = radius;
		this.albedo = albedo;
		this.orbit = orbit;
	}
	public override string ToString(){
		return "Planet {\n\tmass = " + Math.Round(this.mass/Constants.earth.mass, 3).ToString() + 
				" M_earth\n\tradius = " + Math.Round(this.radius/Constants.earth.radius, 2).ToString() + 
				" R_earth\n\talbedo = " + Math.Round(this.albedo, 3).ToString() + 
				"\n}";
	}
}
public class Orbit
{
	public Star star;
	public double aop, ecc, man, sma;
	public Orbit(Star star, double aop, double ecc, double man, double sma){
		this.star = star;
		this.aop = aop;
		this.ecc = ecc;
		this.man = man;
		this.sma = sma;
	}
	public override string ToString(){
		return "Orbit {\n\tsma " + (this.sma / Constants.au).ToString() + " au\n}";
	}
	public double Apoapsis(){
		return (1+this.ecc)*this.sma; 
	}
	public double[] Cartesian(double time){ // confirmed to work 100%
		double E = this.EccentricAnomaly(time);
		double nu = this.TrueAnomaly(time);
		double r_c = this.sma * (1 - this.ecc * Math.Cos(E));
		return new double[] {r_c * Math.Cos(nu), r_c * Math.Sin(nu)};
	}
	public double EccentricAnomaly(double time){ // confirmed to work 100%
		float tol = 1e-10F;
		double M = (this.man + 2*Math.PI*time/this.Period()) % (2*Math.PI);
		double E = M;
		double oldE = E + 2*tol;
		while (tol < Math.Abs(E-oldE)){
			oldE = E;
			E = M + this.ecc * Math.Sin(E);
		}
		return E;
	}
	public double Period(){
		return 2 * Math.PI * Math.Sqrt(Math.Pow(this.sma, 3) / this.star.Mu());
	}
	public double TrueAnomaly(double time){
		double E = this.EccentricAnomaly(time);
		return 2 * Math.Atan2(Math.Sqrt(1+this.ecc) * Math.Sin(E/2), Math.Sqrt(1-this.ecc) * Math.Cos(E/2));
	}
}
public class Star : Body
{
	public double luminosity;
	public ushort temperature;
	public Star(double mass, double radius, double luminosity, ushort temperature) : base(mass, radius){
		this.mass = mass;
		this.radius = radius;
		this.luminosity = luminosity;
		this.temperature = temperature;
	}
	public override string ToString(){
		return "Star {\n\tmass = " + Math.Round(this.mass/Constants.sun.mass, 2).ToString() + 
				" M_sun\n\tradius = " + Math.Round(this.radius/Constants.sun.radius, 2).ToString() + 
				" R_sun\n\tluminosity = " + Math.Round(this.luminosity/Constants.sun.luminosity, 5).ToString() + 
				" L_sun\n\ttemperature = " + this.temperature.ToString() + 
				" K\n}";
	}
	public double FrostLine(){
		return 3*Constants.au * Math.Sqrt(this.luminosity / Constants.sun.luminosity);
	}
}
public class StarSystem
{
	public Star primary;
	public Planet[] secondaries;
	public StarSystem(Star primary, Planet[] secondaries){
		this.primary = primary;
		this.secondaries = secondaries;
	}
	public override string ToString(){
		string o = "System {\n\t" + this.primary.ToString();
		for (byte i=0; i < this.secondaries.Length; i++){
			o += "\n" + this.secondaries[i].ToString();
			o += "\n" + this.secondaries[i].orbit.ToString();
		}
		return o + "\n}";
	}
	public double MaxSMA(){
		return this.secondaries[this.secondaries.Length-1].orbit.Apoapsis();
	}
	public Bitmap Map(ushort size){
		byte orbitResolution = 64;
		Bitmap bitmap = new Bitmap(size, size);
		Graphics g = Graphics.FromImage(bitmap);
		g.Clear(Color.Black);
		// print star!
		Color starColor = Color.Yellow; // new Color(255, 255, 255, 0); // ARGB
		Brush starBrush = new SolidBrush(starColor);
		// bitmap.SetPixel(size/2, size/2, starColor);
		g.FillEllipse(starBrush, size/2-3, size/2-3, 5, 5);
		// compute orbits!
		ushort[,,] pointList = new ushort[this.secondaries.Length, orbitResolution, 2]; // pointlist[planet][time][x/y]
		for (byte i=0; i < this.secondaries.Length; i++){
			double period = this.secondaries[i].orbit.Period();
			for (byte j=0; j < orbitResolution; j++){
				double[] absCoords = this.secondaries[i].orbit.Cartesian(period * j / orbitResolution);
				double absx = absCoords[0];
				double absy = absCoords[1];
				double a = this.MaxSMA();
				pointList[i, j, 0] = (ushort)Constants.Remap(absx, new double[] {-a, a}, new double[] {0, size});
				pointList[i, j, 1] = (ushort)Constants.Remap(absy, new double[] {-a, a}, new double[] {0, size});
			}
		}
		// print orbits!
		Color orbitColor = Color.Red;
		Pen orbitPen = new Pen(orbitColor);
		for (byte i=0; i < this.secondaries.Length; i++){
			for (byte j=0; j < orbitResolution; j++){
				ushort x0 = pointList[i, j, 0];
				ushort y0 = pointList[i, j, 1];
				int new_j = j+1 < orbitResolution ? j+1 : 0;
				ushort x1 = pointList[i, new_j, 0];
				ushort y1 = pointList[i, new_j, 1];
				g.DrawLine(orbitPen, x0, y0, x1, y1);
			}
		}
		// print planets!
		Color planetColor = Color.White;
		Brush planetBrush = new SolidBrush(planetColor);
		for (byte i=0; i < this.secondaries.Length; i++){
			double[] absCoords = this.secondaries[i].orbit.Cartesian(0);
			double absx = absCoords[0];
			double absy = absCoords[1];
			double a = this.MaxSMA();
			ushort x = (ushort)Constants.Remap(absx, new double[] {-a, a}, new double[] {0, size});
			ushort y = (ushort)Constants.Remap(absy, new double[] {-a, a}, new double[] {0, size});
			// bitmap.SetPixel(x, y, planetColor);
			g.FillEllipse(planetBrush, x-2, y-2, 3, 3);
		}
		return bitmap;
	}
}
public class Interface : Form
{
	public Button button1, planetButton, printButton;
	public Label planetSelectorLabel;
	public NumericUpDown planetSelector;
	public PictureBox systemMap;
	public TableLayoutPanel interfaceTable, overTable;
	public Interface(){
		this.Size = new Size(375, 505); // width, height
		this.Text = "Mocha's Stargen";
		this.AutoSize = true;
		// overtable
		overTable = new TableLayoutPanel();
		overTable.Dock = DockStyle.Fill;
		overTable.CellBorderStyle = TableLayoutPanelCellBorderStyle.Single;
		this.Controls.Add(overTable);
		// main table
		interfaceTable = new TableLayoutPanel();
		interfaceTable.Dock = DockStyle.Fill;
		interfaceTable.CellBorderStyle = TableLayoutPanelCellBorderStyle.Single;
		overTable.Controls.Add(interfaceTable, 0, 0);
		// star button
		button1 = new Button();
		button1.Text = "View Star";
		button1.Click += new EventHandler(button1_Click);
		interfaceTable.Controls.Add(button1, 0, 1);
		// planet selector label
		planetSelectorLabel = new Label();
		planetSelectorLabel.Text = "Planet ID:";
		interfaceTable.Controls.Add(planetSelectorLabel, 0, 0);
		// planet selector
		planetSelector = new NumericUpDown();
		planetSelector.Minimum = 0;
		planetSelector.Maximum = Program.system.secondaries.Length - 1;
		interfaceTable.Controls.Add(planetSelector, 1, 0);
		// planet button
		planetButton = new Button();
		planetButton.Text = "View Planet";
		planetButton.Click += new EventHandler(planetButtonClick);
		interfaceTable.Controls.Add(planetButton, 2, 0);
		// print button
		printButton = new Button();
		printButton.Text = "Export System";
		printButton.Width = 100;
		printButton.Click += new EventHandler(printButtonClick);
		interfaceTable.Controls.Add(printButton, 2, 1);
		// map
		systemMap = new PictureBox();
		systemMap.Image = Program.system.Map(350);
		systemMap.Size = new Size(350, 350);
		overTable.Controls.Add(systemMap, 0, 1);
		// interfaceTable.SetColumnSpan(interfaceTable, 3);
		// interfaceTable.Columns[0].Width = 100;
	}
	private void button1_Click(object sender, EventArgs e){
		MessageBox.Show(Program.system.primary.ToString(), "Star");
	}
	private void planetButtonClick(object sender, EventArgs e){
		byte pid = (byte)planetSelector.Value;
		MessageBox.Show(Program.system.secondaries[pid].ToString(), "Planet " + pid.ToString());
	}
	private void printButtonClick(object sender, EventArgs e){
		Program.system.Map(1200).Save(@"export.png", System.Drawing.Imaging.ImageFormat.Png);
		File.WriteAllText("export.txt", Program.system.ToString());
	}
}