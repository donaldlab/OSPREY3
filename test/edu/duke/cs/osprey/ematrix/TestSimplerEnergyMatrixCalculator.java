package edu.duke.cs.osprey.ematrix;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;

import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class TestSimplerEnergyMatrixCalculator extends TestBase {
	
	@Test
	public void discreteALAat10() {
		SimpleConfSpace confSpace = makeConfSpace(false, 10, "ALA");
		assertSingles(confSpace, new double[] {
			-12.69666891725211
		});
	}
	@Test
	public void discreteALAat11() {
		SimpleConfSpace confSpace = makeConfSpace(false, 11, "ALA");
		assertSingles(confSpace, new double[] {
			-13.15430764159240
		});
	}
	@Test
	public void discreteALAat12() {
		SimpleConfSpace confSpace = makeConfSpace(false, 12, "ALA");
		assertSingles(confSpace, new double[] {
			-11.45813233944024
		});
	}
	@Test
	public void discreteALAat13() {
		SimpleConfSpace confSpace = makeConfSpace(false, 13, "ALA");
		assertSingles(confSpace, new double[] {
			-11.87416711976494
		});
	}

	@Test
	public void discreteVALat10() {
		SimpleConfSpace confSpace = makeConfSpace(false, 10, "VAL");
		assertSingles(confSpace, new double[] {
			396.25083612502620,-14.69867910000487,174.46446174980517
		});
	}
	@Test
	public void discreteVALat11() {
		SimpleConfSpace confSpace = makeConfSpace(false, 11, "VAL");
		assertSingles(confSpace, new double[] {
			19.47592126111716,-17.17516256426266,31.48218026243530
		});
	}
	@Test
	public void discreteVALat12() {
		SimpleConfSpace confSpace = makeConfSpace(false, 12, "VAL");
		assertSingles(confSpace, new double[] {
			9.34378533431910,-12.63153233953821,7.96157591033511
		});
	}
	@Test
	public void discreteVALat13() {
		SimpleConfSpace confSpace = makeConfSpace(false, 13, "VAL");
		assertSingles(confSpace, new double[] {
			79.86024836692275,40.47258842584643,36.78125729570677
		});
	}

	@Test
	public void discreteARGat10() {
		SimpleConfSpace confSpace = makeConfSpace(false, 10, "ARG");
		assertSingles(confSpace, new double[] {
			56343048.23570737000000,579.39706538145960,81164969.02657376000000,655.14863804486970,5314096159.72615050000000,297.48056280318497,
			73902464.83018443000000,-19.02521468677143,2050.49278161737450,21244.43943755320300,128376.77398585914000,253061.91563428318000,
			114.29837097221494,221.24780273769570,-35.19546263798888,-35.13999597466904,10799.96670343954800,-28.85857834698063,
			-34.54118356269152,-34.64087003965189,-33.37435942873476,-33.32153284806847,-31.07853882634830,-31.43646074063365,
			-33.63574214824589,-33.62072538908934,482.99340738297460,290.91410026715960,87.21343586034403,922.38341419046650,
			1199.63536823797450,7088.21938326955750,1034.39413561292700,565.77124434439120
		});
	}
	@Test
	public void discreteARGat11() {
		SimpleConfSpace confSpace = makeConfSpace(false, 11, "ARG");
		assertSingles(confSpace, new double[] {
			102715827.07099618000000,4478729.76795961150000,1088314400919.83830000000000,4731987229.52829100000000,118811045.67737170000000,2251109.50034663800000,
			13969474222.74331700000000,711084.18780993320000,1434663451.57182930000000,67224786.93445627000000,36794.39342142169000,50497.72302602212600,
			148494.69699242053000,813757081027.37110000000000,2507.58467397729870,479669753.20372080000000,263473.82004322980000,16490651.01640772400000,
			101413.43250321630000,397687303.71608310000000,559114.83296944340000,337808232.32578045000000,156725.08209436748000,10188.78319922956000,
			19360963590.30088000000000,2842498.36887226860000,43831237523.45173000000000,878705.05301469860000,1919101.85156615820000,1702945271.08319900000000,
			17896373616.69212700000000,394095.15783572575000,7467989.73984597250000,385453411.08065370000000
		});
	}
	@Test
	public void discreteARGat12() {
		SimpleConfSpace confSpace = makeConfSpace(false, 12, "ARG");
		assertSingles(confSpace, new double[] {
			1541248512.46164510000000,141215726.53431967000000,280.69261645517090,265.71960888584374,1178912671.10312180000000,109.49300692249138,
			111.09336466811246,-31.24584640732834,-31.82027384452459,-19.84327380071865,-32.97299894058452,-31.29614571750964,
			-33.96830425277512,107.84524374929426,-33.99084926119196,-34.49114881050276,-25.06817825043814,8835.45014702011000,
			-23.22641409654541,-36.10774173666279,110.15869015804905,110.64623789174462,110.07761966910849,95346.09233818240000,
			2834.84455996298950,821.60620153520340,12786043.78687814600000,23849789.71948394600000,18986106.87954741700000,5679.31416558902600,
			1723.57148204747840,2576.11879220284500,-32.30515069986586,-32.56288830883562
		});
	}
	@Test
	public void discreteARGat13() {
		SimpleConfSpace confSpace = makeConfSpace(false, 13, "ARG");
		assertSingles(confSpace, new double[] {
			3337572.77991629950000,3294152.54535851500000,710093302.51250680000000,15917913.71498998600000,133099709055.01216000000000,4801563729505.82800000000000,
			4801576419413.77800000000000,533038439.79168690000000,443126.72956875403000,226846.41225531689000,20837774.89774106000000,3234924254.02268360000000,
			3068019180.16881900000000,3335653044.78905770000000,413248843.29690343000000,222016.47451258264000,442682.27509270210000,76626.97231332821000,
			133323.68481809003000,882227184.28396300000000,56003204.98683555000000,35825720.59591460000000,2180923.73791479830000,499476.55813721340000,
			702999.50943198740000,6653.66331217276100,97957.17366154466000,322028.69703886006000,131.93949899186290,60142219646.89902500000000,
			1874816.58250321980000,15988.04109407721600,73635143.88237545000000,159471081260320.72000000000000
		});
	}

	@Test
	public void discretePROat10() {
		SimpleConfSpace confSpace = makeConfSpace(false, 10, "PRO");
		assertSingles(confSpace, new double[] {
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY
		});
	}
	@Test
	public void discretePROat11() {
		SimpleConfSpace confSpace = makeConfSpace(false, 11, "PRO");
		assertSingles(confSpace, new double[] {
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY
		});
	}
	@Test
	public void discretePROat12() {
		SimpleConfSpace confSpace = makeConfSpace(false, 12, "PRO");
		assertSingles(confSpace, new double[] {
			101771.19435558538000,17766.72158716243000
		});
	}
	@Test
	public void discretePROat13() {
		SimpleConfSpace confSpace = makeConfSpace(false, 13, "PRO");
		assertSingles(confSpace, new double[] {
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY
		});
	}

	@Test
	public void discreteLEUat10() {
		SimpleConfSpace confSpace = makeConfSpace(false, 10, "LEU");
		assertSingles(confSpace, new double[] {
			4694555.93770604400000,17.48761152046505,37.97563037855740,-12.68443340951747,1882.42275026271820
		});
	}
	@Test
	public void discreteLEUat11() {
		SimpleConfSpace confSpace = makeConfSpace(false, 11, "LEU");
		assertSingles(confSpace, new double[] {
			58784.21884701415500,51.80506084516539,60.34644321090021,407.64872336278670,1096.17459418563300
		});
	}
	@Test
	public void discreteLEUat12() {
		SimpleConfSpace confSpace = makeConfSpace(false, 12, "LEU");
		assertSingles(confSpace, new double[] {
			3428844.76079163750000,-12.23363676601438,-1.94154544513084,29.20342022336552,379.54234151032900
		});
	}
	@Test
	public void discreteLEUat13() {
		SimpleConfSpace confSpace = makeConfSpace(false, 13, "LEU");
		assertSingles(confSpace, new double[] {
			20666744.72983168800000,1199604.68333470900000,2131302.33510759660000,15.91099873397199,9.32540533793335
		});
	}

	@Test
	public void discretePHEat10() {
		SimpleConfSpace confSpace = makeConfSpace(false, 10, "PHE");
		assertSingles(confSpace, new double[] {
			162516655.03290004000000,1161.22928139693500,187.95665040172767,2623.45743171848650
		});
	}
	@Test
	public void discretePHEat11() {
		SimpleConfSpace confSpace = makeConfSpace(false, 11, "PHE");
		assertSingles(confSpace, new double[] {
			305984828.89335220000000,13652.82094287018700,2972869.03802562130000,4545.93569607702600
		});
	}
	@Test
	public void discretePHEat12() {
		SimpleConfSpace confSpace = makeConfSpace(false, 12, "PHE");
		assertSingles(confSpace, new double[] {
			5917472.55506810450000,30.60868320878065,2018.30666331423480,1009468.03169236320000
		});
	}
	@Test
	public void discretePHEat13() {
		SimpleConfSpace confSpace = makeConfSpace(false, 13, "PHE");
		assertSingles(confSpace, new double[] {
			1288151813495.98830000000000,793844506.96673920000000,17792.36400344047000,61.11637413764255
		});
	}

	@Test
	public void discreteTHRat10() {
		SimpleConfSpace confSpace = makeConfSpace(false, 10, "THR");
		assertSingles(confSpace, new double[] {
			-9.60831835704874,-9.80499482321600,-8.44916624055827,-9.53404243956254,-8.83633582675501,-8.92348661381635,
			477.58605421839360,477.04302629278584,477.60889815208077,477.28505954272356,477.32202831556700,477.67549535206580,
			-18.21692687254067,-18.21235738221050,-17.70691838446186,-18.29755726431351,-17.83234657656701,-18.06465403147951
		});
	}
	@Test
	public void discreteTHRat11() {
		SimpleConfSpace confSpace = makeConfSpace(false, 11, "THR");
		assertSingles(confSpace, new double[] {
			-16.29018735981431,-17.17176568515227,-16.52469102070694,-16.96101436915492,-16.80946785201232,-16.23786712965228,
			15.80240534128683,15.85635687353529,15.68550479881661,15.99313207113857,15.79783477886055,15.57354960211940,
			-19.16168199617284,-19.47303724751162,-19.36171080933660,-19.33421558730525,-19.47446926042962,-19.28299834547580
		});
	}
	@Test
	public void discreteTHRat12() {
		SimpleConfSpace confSpace = makeConfSpace(false, 12, "THR");
		assertSingles(confSpace, new double[] {
			-13.44999350724677,-14.14201918723017,-13.41215613829489,-13.95450332480944,-13.31846689403499,-13.56523767004068,
			7.06645091052562,6.74514352271787,6.45605813261890,7.10090816568600,6.50411226061839,6.68839660073257,
			-15.50252637106665,-15.74943617747870,-16.77782634824157,-15.47737563471156,-16.46895516436925,-15.98734340905233
		});
	}
	@Test
	public void discreteTHRat13() {
		SimpleConfSpace confSpace = makeConfSpace(false, 13, "THR");
		assertSingles(confSpace, new double[] {
			-11.20778602586725,-10.99698649688516,-10.57845348813976,-10.58078759713326,-10.75289184203710,-11.00224978205113,
			69.95418627108724,71.22615811986262,70.06846973272570,71.08173040894454,70.56826708684080,69.43469194735525,
			38.70681941817195,37.58750579848925,38.26118190432410,38.16983254847784,37.80956147961889,38.55764814467396
		});
	}

	@Test
	public void discreteLYSat10() {
		SimpleConfSpace confSpace = makeConfSpace(false, 10, "LYS");
		assertSingles(confSpace, new double[] {
			676189.83497802600000,15700208.41864713700000,945.76536615675160,563056.44977351280000,555.63885675010950,54.00158282017770,
			123.12015479338420,4.55318023616803,7686.90531851778700,12.31009975227060,-20.30202976208145,-20.10722368461198,
			-14.00727215121759,-19.75679275486062,-19.50929125513042,185.07373878291048,-18.11434795387376,-18.35234067875615,
			-17.61447452455994,-18.60838319454741,-18.19676017757567,234.98723627452216,3334.34941791264650,1337.39417498655700,
			1652.86545385348450,1476.28611823082430,615.59329546531510
		});
	}
	@Test
	public void discreteLYSat11() {
		SimpleConfSpace confSpace = makeConfSpace(false, 11, "LYS");
		assertSingles(confSpace, new double[] {
			2144965.57933704930000,129183007.33876583000000,70908.28908200457000,45059.55821341886600,53160.24495681152000,49221.78169653419000,
			11460.99515413466700,3060.00097405383030,8983.64232341394600,606851.74020416280000,3822.88605749027650,171128349.10491928000000,
			1780.40443988446440,490.64157050183660,776.68693936090160,4072122.11817667630000,1537.36267085065220,895.67628367805530,
			1657.41669715943340,500373278.66949135000000,20094.98729142180700,81836.36155140829000,629712.32464023820000,12800.93160822744400,
			1032824.45681229840000,23826.77465878245000,9344.13717786387000
		});
	}
	@Test
	public void discreteLYSat12() {
		SimpleConfSpace confSpace = makeConfSpace(false, 12, "LYS");
		assertSingles(confSpace, new double[] {
			5174818.03025382700000,930.59627364274250,711.41973654480240,1655997.14317585710000,118.88366757164212,3.14105663448442,
			-18.31914058918627,-17.69042237112890,-17.33883667465622,-19.68179780023592,-19.03361404171609,-19.92052605861886,
			6.33590828855181,-20.96033603806149,-18.93012481539868,559873.65588798290000,190.08833695787650,187.95551590606618,
			654085.89559014560000,5020.92947056868250,933.65631818471150,11938.80878919721500,178446.29320344923000,-16.34979255604794,
			-16.39555591632479,326.25331974635156,-16.57051160317388
		});
	}
	@Test
	public void discreteLYSat13() {
		SimpleConfSpace confSpace = makeConfSpace(false, 13, "LYS");
		assertSingles(confSpace, new double[] {
			2495073.31923758750000,32572019.83250775000000,32414034.73555890000000,31988884.98165150000000,671593843665667.00000000000000,2662956.67390299000000,
			13044118.98643033200000,837151.46416285000000,19740596267.75419600000000,13718260368413210.00000000000000,20134416.50057526000000,160444.17239661064000,
			155849.51713175443000,92822.04300561998000,417019.24560419600000,32593321.72972560700000,531545.51950390020000,11458.16385466861900,
			262.71405188482660,223.63463986832272,107.97807463512346,69.76512955010027,-21.62061547667849,183.31514105711267,
			4271.38763359589800,3642.85786653818900,9369523.44418660400000
		});
	}

	@Test
	public void discreteGLUat10() {
		SimpleConfSpace confSpace = makeConfSpace(false, 10, "GLU");
		assertSingles(confSpace, new double[] {
			325.44708478173890,202183.46753074962000,-5.73448895655529,-25.02723262822549,55553.12585336772000,154.14306453148822,
			-24.85080716659156,14022.48723115565700
		});
	}
	@Test
	public void discreteGLUat11() {
		SimpleConfSpace confSpace = makeConfSpace(false, 11, "GLU");
		assertSingles(confSpace, new double[] {
			26894.25457700122600,483.84721405373625,32.86815391331179,141.88204287536877,193226.33748605152000,20140.29101701810600,
			852.16491572469300,849.72752306572530
		});
	}
	@Test
	public void discreteGLUat12() {
		SimpleConfSpace confSpace = makeConfSpace(false, 12, "GLU");
		assertSingles(confSpace, new double[] {
			187.23749389993240,2337475.66478547870000,-18.86254488153889,-19.22560221544153,9573.98639943477400,2737.92972976652200,
			365.17150205459050,-21.36629730247750
		});
	}
	@Test
	public void discreteGLUat13() {
		SimpleConfSpace confSpace = makeConfSpace(false, 13, "GLU");
		assertSingles(confSpace, new double[] {
			5243.58334916048350,2098.44482079887500,2716.38075705549200,1900.27621655323220,545800.25783020510000,2183.05152739109460,
			-19.13572754956669,124.95719792719579
		});
	}

	@Test
	public void continuousALAat10() {
		SimpleConfSpace confSpace = makeConfSpace(true, 10, "ALA");
		assertSingles(confSpace, new double[] {
			-12.69666891725211
		});
	}
	@Test
	public void continuousALAat11() {
		SimpleConfSpace confSpace = makeConfSpace(true, 11, "ALA");
		assertSingles(confSpace, new double[] {
			-13.15430764159240
		});
	}
	@Test
	public void continuousALAat12() {
		SimpleConfSpace confSpace = makeConfSpace(true, 12, "ALA");
		assertSingles(confSpace, new double[] {
			-11.45813233944024
		});
	}
	@Test
	public void continuousALAat13() {
		SimpleConfSpace confSpace = makeConfSpace(true, 13, "ALA");
		assertSingles(confSpace, new double[] {
			-11.87416711976494
		});
	}

	@Test
	public void continuousVALat10() {
		SimpleConfSpace confSpace = makeConfSpace(true, 10, "VAL");
		assertSingles(confSpace, new double[] {
			62.14496323265180,-14.90025408831916,32.32795959778965
		});
	}
	@Test
	public void continuousVALat11() {
		SimpleConfSpace confSpace = makeConfSpace(true, 11, "VAL");
		assertSingles(confSpace, new double[] {
			12.63690402989719,-18.75742019128871,13.23250240028139
		});
	}
	@Test
	public void continuousVALat12() {
		SimpleConfSpace confSpace = makeConfSpace(true, 12, "VAL");
		assertSingles(confSpace, new double[] {
			8.02165581418594,-12.65214486629394,7.92548474344134
		});
	}
	@Test
	public void continuousVALat13() {
		SimpleConfSpace confSpace = makeConfSpace(true, 13, "VAL");
		assertSingles(confSpace, new double[] {
			65.84547279228687,5.31719889505351,23.87885411379908
		});
	}

	@Test
	public void continuousARGat10() {
		SimpleConfSpace confSpace = makeConfSpace(true, 10, "ARG");
		assertSingles(confSpace, new double[] {
			16498.85853555735000,49.61991066715400,1225786.35400247080000,275.47092378955550,115008.04656973168000,73.98644576633401,
			284793.47389117075000,-33.99270338766021,-31.63326797832978,3.51840546461714,70.83866137230139,323.80170144963460,
			-33.54615326415861,-32.02649044870589,-35.72394550350501,-35.15700003650038,-33.00617441735968,-33.49645383264431,
			-34.83627895759657,-35.34717155809033,-34.58383593723006,-34.29028002165823,-33.21714053835623,-34.10877582285354,
			-34.50232546824618,-34.80507666982112,-27.49539955822932,-26.84787866879763,-31.47581464840185,12.76493870737142,
			16.37041711339186,154.05995282636067,2.10379272394659,-26.85734991012119
		});
	}
	@Test
	public void continuousARGat11() {
		SimpleConfSpace confSpace = makeConfSpace(true, 11, "ARG");
		assertSingles(confSpace, new double[] {
			418338.18846395790000,41610.85452582272000,170468.25651668722000,109114.59593332998000,101595.06165900196000,6560.91115104822000,
			280669.35573200276000,3586.31427188826600,6667.44944008919600,28555.72545483618500,1947.74360517813650,1145.02563194562780,
			4460.22851833337400,81360.12228444825000,475.46461644986330,37527169.24045370000000,648.92812727092700,4716.23175377909300,
			67837.32525598740000,102065.41524319783000,3553.78461265236700,7029.74926664652100,9174.35497460325300,7043.22437057653400,
			4824.75086131642100,972.26307063767070,208242.01387457742000,28022.48318397822700,29256.23911532655600,6243.39438827056300,
			288685.08257663855000,5040.88017666709400,3510.69180320959140,19023.85343407539000
		});
	}
	@Test
	public void continuousARGat12() {
		SimpleConfSpace confSpace = makeConfSpace(true, 12, "ARG");
		assertSingles(confSpace, new double[] {
			62148421.25385725000000,7103.96778424942200,27.98290313257352,25.83709695620482,3993475.73945789640000,8.53400706942282,
			10.23931965446292,-32.97205666364494,-33.41163983503019,-34.87302262858511,-33.57385399337947,-33.92747973291882,
			-34.85448646280618,-28.31912487195003,-34.31051784324831,-35.03795640058898,-33.97111575091919,402.93025905703920,
			-34.42472938998923,-37.08260533574643,-17.50319551499512,-20.04789932530111,-17.72834533697541,2003.43180795270060,
			14.57040977786101,9.88212136637196,291.67532788886166,463.52910780965740,210239.48788920890000,-31.32759144692056,
			-31.07466353165080,-20.14991812915281,-33.45231196774338,-33.55167955177462
		});
	}
	@Test
	public void continuousARGat13() {
		SimpleConfSpace confSpace = makeConfSpace(true, 13, "ARG");
		assertSingles(confSpace, new double[] {
			5625.00036988318900,1351.66219006210080,23185219.53341683000000,150081.45108063697000,3725698.44569358000000,1397059.84851077200000,
			2360803.89811731830000,1256768.92289058030000,5290.94042040810000,28587.68018937310300,1000260.12100091860000,407594.48337878570000,
			576470.15223369190000,127582.03425339272000,764945.04828077350000,12213.44282904487200,47083.96406633398000,377.59792096107634,
			20231.78855272149300,7341984.26475459700000,147320.39693789775000,515612.30476869230000,747.43606611127830,29438.71288985778300,
			27976.89114376512600,1180.27380974506830,115.52562061751695,736.35002252199260,-15.56318959203869,71629.93155518571000,
			1985.67403736915620,516.28916574134000,205724.63643789222000,6157162.79574698000000
		});
	}

	@Test
	public void continuousPROat10() {
		SimpleConfSpace confSpace = makeConfSpace(true, 10, "PRO");
		assertSingles(confSpace, new double[] {
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY
		});
	}
	@Test
	public void continuousPROat11() {
		SimpleConfSpace confSpace = makeConfSpace(true, 11, "PRO");
		assertSingles(confSpace, new double[] {
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY
		});
	}
	@Test
	public void continuousPROat12() {
		// NOTE: this is the only test where proline doesn't have conf problems
		// ie, we get finite energies here
		SimpleConfSpace confSpace = makeConfSpace(true, 12, "PRO");
		assertSingles(confSpace, new double[] {
			101771.19435558538000,17766.72158716243000
		});
	}
	@Test
	public void continuousPROat13() {
		SimpleConfSpace confSpace = makeConfSpace(true, 13, "PRO");
		assertSingles(confSpace, new double[] {
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY
		});
	}

	@Test
	public void continuousLEUat10() {
		SimpleConfSpace confSpace = makeConfSpace(true, 10, "LEU");
		assertSingles(confSpace, new double[] {
			837097.50184635750000,-12.64708073819766,-11.61589985913956,-14.32694917861948,31.30244943565061
		});
	}
	@Test
	public void continuousLEUat11() {
		SimpleConfSpace confSpace = makeConfSpace(true, 11, "LEU");
		assertSingles(confSpace, new double[] {
			18135.04661794376000,37.42861587418536,-10.38732315034924,165.34820435572357,617.96253281987050
		});
	}
	@Test
	public void continuousLEUat12() {
		SimpleConfSpace confSpace = makeConfSpace(true, 12, "LEU");
		assertSingles(confSpace, new double[] {
			72769.43689368735000,-13.39863599385466,-11.81158022790220,-8.22535905852000,60.62288888313260
		});
	}
	@Test
	public void continuousLEUat13() {
		SimpleConfSpace confSpace = makeConfSpace(true, 13, "LEU");
		assertSingles(confSpace, new double[] {
			29597.73285663485500,11986.15872043812300,68168.69646346256000,-0.25051398177132,-13.95447870698099
		});
	}

	@Test
	public void continuousPHEat10() {
		SimpleConfSpace confSpace = makeConfSpace(true, 10, "PHE");
		assertSingles(confSpace, new double[] {
			3637582.30473935980000,-2.16394976228906,-12.96626948159071,64.35642712008150
		});
	}
	@Test
	public void continuousPHEat11() {
		SimpleConfSpace confSpace = makeConfSpace(true, 11, "PHE");
		assertSingles(confSpace, new double[] {
			1084876.06144880460000,10546.71879268759600,8144.65059196565000,969.91257233874660
		});
	}
	@Test
	public void continuousPHEat12() {
		SimpleConfSpace confSpace = makeConfSpace(true, 12, "PHE");
		assertSingles(confSpace, new double[] {
			584646.64827945140000,-8.08549614187687,458.28353876499290,7539.20871500066000
		});
	}
	@Test
	public void continuousPHEat13() {
		SimpleConfSpace confSpace = makeConfSpace(true, 13, "PHE");
		assertSingles(confSpace, new double[] {
			953367.40082424300000,2669428.07767566200000,529.41028382188030,45.16556412329300
		});
	}

	@Test
	public void continuousTHRat10() {
		SimpleConfSpace confSpace = makeConfSpace(true, 10, "THR");
		assertSingles(confSpace, new double[] {
			-14.14088218030851,-13.61337500568927,-12.47648721137899,-13.47526078237647,-12.96280392316086,-13.66364404131135,
			72.69589104943117,72.16024188773366,72.66725077657809,72.42459278241694,72.35430562544755,72.78646167240302,
			-18.63082757062709,-18.66913014953231,-18.10865163518616,-18.74689911925026,-18.34362665140533,-18.48544204570376
		});
	}
	@Test
	public void continuousTHRat11() {
		SimpleConfSpace confSpace = makeConfSpace(true, 11, "THR");
		assertSingles(confSpace, new double[] {
			-16.41119821919324,-17.38433438338003,-16.65722556618435,-17.35006906099130,-16.90102105634410,-16.33061311605634,
			11.17087226794999,11.25696932568680,11.08371954543446,11.35839150735534,11.19780109631335,11.00637118313378,
			-20.35670844313994,-20.69418442738629,-20.53221585009801,-20.56563949846909,-20.67147501721451,-20.44466223129886
		});
	}
	@Test
	public void continuousTHRat12() {
		SimpleConfSpace confSpace = makeConfSpace(true, 12, "THR");
		assertSingles(confSpace, new double[] {
			-13.67434642734197,-14.36620876030489,-13.61721889566467,-14.17435405084187,-13.48060293414165,-13.79109621707043,
			4.75823761475827,4.49854224135211,4.22479170987210,4.87135721379240,4.26631250636086,4.41490470307813,
			-15.94888610352922,-16.21948455764298,-17.29081966249179,-15.91307017871355,-16.93553543549995,-16.56807799714450
		});
	}
	@Test
	public void continuousTHRat13() {
		SimpleConfSpace confSpace = makeConfSpace(true, 13, "THR");
		assertSingles(confSpace, new double[] {
			-15.51592359322734,-14.96467833787075,-14.92796981106832,-15.22631931808176,-14.99044249209065,-15.37336894207380,
			28.08513250719079,29.49013599908995,28.27151763761004,29.14973503101176,28.91176772988028,27.61144827026070,
			0.87955141918212,-0.36122928884290,0.53333311778052,0.25279936993078,-0.02741562202855,0.83639591457118
		});
	}

	@Test
	public void continuousLYSat10() {
		SimpleConfSpace confSpace = makeConfSpace(true, 10, "LYS");
		assertSingles(confSpace, new double[] {
			122.33018659764190,8829.28402244064000,639.12268064468890,18249.60813105754300,185.31027899969908,-10.77066497786510,
			-14.52306296387104,-19.19193720496622,42.58529548785467,-18.03329368103528,-20.67996132522209,-20.41799051727656,
			-19.43400042197833,-20.34966374019472,-20.38572483883375,-18.40248602870868,-19.62646267505800,-19.65180850357086,
			-19.26049660241277,-19.89748992258656,-19.69114453964284,-9.31454217997359,3.55900001611279,35.47312474979685,
			42.97601282125492,176.94734402042230,-11.36474954589414
		});
	}
	@Test
	public void continuousLYSat11() {
		SimpleConfSpace confSpace = makeConfSpace(true, 11, "LYS");
		assertSingles(confSpace, new double[] {
			40000.67635194762000,849300.77472173180000,21240.86608054120400,14070.37462520016600,10515.13862040789400,484.77645457687550,
			124.38506909971699,61.89267324284393,1789.09848622818440,959.18937936563900,582.70961188189000,10741723.29623147000000,
			650.36294381939920,150.20964846105713,62.03454799944217,1247.29865813433500,82.06437822649137,396.57529014313290,
			867.06571400385420,8678.48374305250000,414.29781428006817,5550.92729156558300,4237.66288424078300,224.76803469136547,
			2114.52053215687300,6456.79755901741900,582.33667174877610
		});
	}
	@Test
	public void continuousLYSat12() {
		SimpleConfSpace confSpace = makeConfSpace(true, 12, "LYS");
		assertSingles(confSpace, new double[] {
			1401.23256798114430,58.83776724569968,53.76554865279540,57942.62377191142000,25.22569498314587,-18.78906508526520,
			-19.02020901528003,-18.93581110673657,-18.53320686547477,-20.97641021084721,-19.61867194608093,-20.30048152740476,
			-18.50426037191262,-21.76221162963867,-19.75034609498356,468.40540576498324,2.95248081880122,1.33935979407349,
			3428.08334064458900,39.48907637933474,36.40375043446331,756.15159588607420,12088.63720169499700,-20.13203501184136,
			-18.46105966689958,-16.41905198854522,-18.36239818354909
		});
	}
	@Test
	public void continuousLYSat13() {
		SimpleConfSpace confSpace = makeConfSpace(true, 13, "LYS");
		assertSingles(confSpace, new double[] {
			1052.77054722525210,175466.92382837873000,117859.57315758846000,355371.12072144094000,71638667.35572597000000,331548.39231622266000,
			85860.96694414206000,26073.93938082118600,9635741.87472551700000,6335745.52947419000000,20402.69552537183700,13554.12279986917700,
			13998.22953026067800,2260.43119275999920,1056.95336999283930,975.35413841175670,247.00766006783735,5696.54111373476500,
			-11.17754725030058,18.60973834870556,-8.67147754035173,-19.19745400975485,-23.33997281599913,22.54423808661629,
			64.20346642219633,93.62695329490340,5446.71957921134800
		});
	}

	@Test
	public void continuousGLUat10() {
		SimpleConfSpace confSpace = makeConfSpace(true, 10, "GLU");
		assertSingles(confSpace, new double[] {
			85.14390853092256,39439.87458391501000,-18.13597559307399,-25.27202980944163,67.61016355691078,-5.34081109009057,
			-25.31741272474850,156.42539750875530
		});
	}
	@Test
	public void continuousGLUat11() {
		SimpleConfSpace confSpace = makeConfSpace(true, 11, "GLU");
		assertSingles(confSpace, new double[] {
			585.86671180958330,183.22538824663758,11.95523111386433,26.96357357281653,1896.96946036113560,1360.78263445343700,
			22.19229815935954,498.57018575353140
		});
	}
	@Test
	public void continuousGLUat12() {
		SimpleConfSpace confSpace = makeConfSpace(true, 12, "GLU");
		assertSingles(confSpace, new double[] {
			20.16380604035198,14374.08858031909400,-20.13068447164413,-19.93315573238242,121.32233932029351,380.87396495108820,
			-5.03900571051331,-21.72184735285783
		});
	}
	@Test
	public void continuousGLUat13() {
		SimpleConfSpace confSpace = makeConfSpace(true, 13, "GLU");
		assertSingles(confSpace, new double[] {
			93.81720795230409,1691.02918971144600,677.57673143633870,1736.49666416634070,5780.45863350600100,304.12738407685400,
			-23.26274146400608,20.26498689961222
		});
	}
	
	@Test
	public void discreteALAtoVAL() {
		SimpleConfSpace confSpace = makeConfSpace(false, "ALA", "THR", "ARG", "VAL");
		assertEnergyMatrix(
			confSpace,
			makeExpectedEmatDiscreteALAtoVAL(confSpace),
			makeEmatCalc(confSpace).calcEnergyMatrix()
		);
	}
	
	@Test
	public void discreteGLYtoGLU() {
		SimpleConfSpace confSpace = makeConfSpace(false, "GLY", "SER", "ASN", "GLU");
		assertEnergyMatrix(
			confSpace,
			makeExpectedEmatDiscreteGLYtoGLU(confSpace),
			makeEmatCalc(confSpace).calcEnergyMatrix()
		);
	}
	
	@Test
	public void discreteASPtoPHE() {
		SimpleConfSpace confSpace = makeConfSpace(false, "ASP", "CYS", "MET", "PHE");
		assertEnergyMatrix(
			confSpace,
			makeExpectedEmatDiscreteASPtoPHE(confSpace),
			makeEmatCalc(confSpace).calcEnergyMatrix()
		);
	}
	
	@Test
	public void discreteGLNtoTRP() {
		SimpleConfSpace confSpace = makeConfSpace(false, "GLN", "ILE", "LEU", "TRP");
		assertEnergyMatrix(
			confSpace,
			makeExpectedEmatDiscreteGLNtoTRP(confSpace),
			makeEmatCalc(confSpace).calcEnergyMatrix()
		);
	}
	
	@Test
	public void discreteHIStoTYR() {
		SimpleConfSpace confSpace = makeConfSpace(false, "HIS", "LYS", "TYR");
		assertEnergyMatrix(
			confSpace,
			makeExpectedEmatDiscreteHIStoTYR(confSpace),
			makeEmatCalc(confSpace).calcEnergyMatrix()
		);
	}
	
	@Test
	public void discretePROtoPRO() {
		SimpleConfSpace confSpace = makeConfSpace(false, "PRO");
		assertEnergyMatrix(
			confSpace,
			makeExpectedEmatDiscretePROtoPRO(confSpace),
			makeEmatCalc(confSpace).calcEnergyMatrix()
		);
	}
	
	@Test
	public void discreteLYStoVAL() {
		SimpleConfSpace confSpace = makeConfSpace(false, "LYS", "PRO", "VAL");
		assertEnergyMatrix(
			confSpace,
			makeExpectedEmatDiscreteLYStoVAL(confSpace),
			makeEmatCalc(confSpace).calcEnergyMatrix()
		);
	}
	
	@Test
	public void continuousALAtoVAL() {
		SimpleConfSpace confSpace = makeConfSpace(true, "ALA", "THR", "ARG", "VAL");
		assertEnergyMatrix(
			confSpace,
			makeExpectedEmatContinuousALAtoVAL(confSpace),
			makeEmatCalc(confSpace).calcEnergyMatrix()
		);
	}

	@Test
	public void continuousGLYtoGLU() {
		SimpleConfSpace confSpace = makeConfSpace(true, "GLY", "SER", "ASN", "GLU");
		assertEnergyMatrix(
			confSpace,
			makeExpectedEmatContinuousGLYtoGLU(confSpace),
			makeEmatCalc(confSpace).calcEnergyMatrix()
		);
	}
	
	@Test
	public void continuousASPtoPHE() {
		SimpleConfSpace confSpace = makeConfSpace(true, "ASP", "CYS", "MET", "PHE");
		assertEnergyMatrix(
			confSpace,
			makeExpectedEmatContinuousASPtoPHE(confSpace),
			makeEmatCalc(confSpace).calcEnergyMatrix()
		);
	}
	
	@Test
	public void continuousGLNtoTRP() {
		SimpleConfSpace confSpace = makeConfSpace(true, "GLN", "ILE", "LEU", "TRP");
		assertEnergyMatrix(
			confSpace,
			makeExpectedEmatContinuousGLNtoTRP(confSpace),
			makeEmatCalc(confSpace).calcEnergyMatrix()
		);
	}
	
	@Test
	public void continuousHIStoTYR() {
		SimpleConfSpace confSpace = makeConfSpace(true, "HIS", "LYS", "TYR");
		assertEnergyMatrix(
			confSpace,
			makeExpectedEmatContinuousHIStoTYR(confSpace),
			makeEmatCalc(confSpace).calcEnergyMatrix()
		);
	}
	
	@Test
	public void continuousPROtoPRO() {
		SimpleConfSpace confSpace = makeConfSpace(true, "PRO");
		assertEnergyMatrix(
			confSpace,
			makeExpectedEmatContinuousPROtoPRO(confSpace),
			makeEmatCalc(confSpace).calcEnergyMatrix()
		);
	}
	
	@Test
	public void continuousLYStoVAL() {
		SimpleConfSpace confSpace = makeConfSpace(true, "LYS", "PRO", "VAL");
		assertEnergyMatrix(
			confSpace,
			makeExpectedEmatContinuousLYStoVAL(confSpace),
			makeEmatCalc(confSpace).calcEnergyMatrix()
		);
	}
	
	private SimpleConfSpace makeConfSpace(boolean doMinimize, String ... aminoAcids) {
		return makeConfSpace(doMinimize, 10, aminoAcids);
	}
	
	private SimpleConfSpace makeConfSpace(boolean doMinimize, int firstResNum, String ... aminoAcids) {
		
		Molecule mol = PDBIO.readFile("examples/python.GMEC/1CC8.ss.pdb");
		Strand strand = new Strand.Builder(mol)
			// explicitly choose Lovell rotamers
			.setTemplateLibrary(new ResidueTemplateLibrary.Builder()
					.clearRotamers()
					.addLovellRotamers()
					.build())
			.build();
		
		int resNum = firstResNum;
		for (String aminoAcid : aminoAcids) {
			strand.flexibility.get("A" + resNum).setLibraryRotamers(aminoAcid);
			if (doMinimize) {
				strand.flexibility.get("A" + resNum).setContinuous();
			}
			resNum++;
		}
		
		return new SimpleConfSpace.Builder().addStrand(strand).build();
	}
	
	private SimplerEnergyMatrixCalculator makeEmatCalc(SimpleConfSpace confSpace) {
		EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setType(EnergyCalculator.Type.CpuOriginalCCD) // use original CCD implementation to match old code energies
			.build();
		return new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
			.build();
	}
	
	private void assertSingles(SimpleConfSpace confSpace, double[] exp) {
		for (int i=0; i<confSpace.positions.get(0).resConfs.size(); i++) {
			assertEnergy(
				confSpace.isContinuouslyFlexible(new RCTuple(0, i)),
				exp[i],
				makeEmatCalc(confSpace).confEcalc.calcSingleEnergy(0, i).energy,
				describe(confSpace, 0, i)
			); 
		}
	}
	
	private static void assertEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix exp, EnergyMatrix obs) {
		
		assertThat(obs.getNumPos(), is(exp.getNumPos()));
		
		for (int pos1=0; pos1<exp.getNumPos(); pos1++) {
			for (int rc1=0; rc1<exp.getNumConfAtPos(pos1); rc1++) {
				
				assertThat(obs.getNumConfAtPos(pos1), is(exp.getNumConfAtPos(pos1)));
				assertSingle(confSpace, exp, obs, pos1, rc1);
				
				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<exp.getNumConfAtPos(pos2); rc2++) {
						
						assertPair(confSpace, exp, obs, pos1, rc1, pos2, rc2);
					}
				}
			}
		}
	}
	
	private static String describe(SimpleConfSpace confSpace, int pos1, int rc1) {
		return String.format("%s %s, %d",
			confSpace.positions.get(pos1).resNum,
			confSpace.positions.get(pos1).resConfs.get(rc1).template.name,
			rc1
		);
	}
	
	private static String describe(SimpleConfSpace confSpace, int pos1, int rc1, int pos2, int rc2) {
		return String.format("%s %s, %d x %s %s, %d",
			confSpace.positions.get(pos1).resNum,
			confSpace.positions.get(pos1).resConfs.get(rc1).template.name,
			rc1,
			confSpace.positions.get(pos2).resNum,
			confSpace.positions.get(pos2).resConfs.get(rc2).template.name,
			rc2
		);
	}

	private static void assertSingle(SimpleConfSpace confSpace, EnergyMatrix exp, EnergyMatrix obs, int pos1, int rc1) {
		assertEnergy(
                        confSpace.isContinuouslyFlexible(new RCTuple(pos1, rc1)),
			exp.getOneBody(pos1, rc1),
			obs.getOneBody(pos1, rc1),
			describe(confSpace, pos1, rc1)
		);
	}
	
	private static void assertPair(SimpleConfSpace confSpace, EnergyMatrix exp, EnergyMatrix obs, int pos1, int rc1, int pos2, int rc2) {
		assertEnergy(
                        confSpace.isContinuouslyFlexible(new RCTuple(pos1, rc1, pos2, rc2)),
			exp.getPairwise(pos1, rc1, pos2, rc2),
			obs.getPairwise(pos1, rc1, pos2, rc2),
			describe(confSpace, pos1, rc1, pos2, rc2)
		);
	}
	
	private static void assertEnergy(boolean isMinimized, double exp, double obs, String description) {
		
		// explicitly check for infinity
		if (exp == Double.POSITIVE_INFINITY) {
			assertThat(description, obs, is(Double.POSITIVE_INFINITY));
		} else {
			
			// for minimized energies, lower observed energy is ok,
			// since improvements to minimizers over time could give us lower energies
			if (isMinimized && obs < exp) {
				return;
			}
			
			// high energies won't be as accurate,
			// so use a looser relative check
			if (Math.abs(exp) > 100) {
				assertThat(description, obs, isRelatively(exp, 1e-4));
				
			// otherwise, use a tight absolute check
			} else {
				assertThat(description, obs, isAbsolutely(exp, 1e-5));
			}
		}
	}
	
	
	// here's some tools to auto-generate the expected energies for this test
	
	public static void main(String[] args) {
		
		TestBase.initDefaultEnvironment();
		
		StringBuilder buf = new StringBuilder();
		for (boolean doMinimize : Arrays.asList(false, true)) {
			
			getExpectedEnergies(buf, doMinimize, "ALA", "THR", "ARG", "VAL");
			getExpectedEnergies(buf, doMinimize, "GLY", "SER", "ASN", "GLU");
			getExpectedEnergies(buf, doMinimize, "ASP", "CYS", "MET", "PHE");
			getExpectedEnergies(buf, doMinimize, "GLN", "ILE", "LEU", "TRP");
			getExpectedEnergies(buf, doMinimize, "HIS", "LYS", "TYR");
			getExpectedEnergies(buf, doMinimize, "PRO");
			getExpectedEnergies(buf, doMinimize, "LYS", "PRO", "VAL");
			
			getExpectedEnergy(buf, doMinimize, "ALA");
			getExpectedEnergy(buf, doMinimize, "VAL");
			getExpectedEnergy(buf, doMinimize, "ARG");
			getExpectedEnergy(buf, doMinimize, "PRO");
			getExpectedEnergy(buf, doMinimize, "LEU");
			getExpectedEnergy(buf, doMinimize, "PHE");
			getExpectedEnergy(buf, doMinimize, "THR");
			getExpectedEnergy(buf, doMinimize, "LYS");
			getExpectedEnergy(buf, doMinimize, "GLU");
		}
		System.out.println(buf);
	}
	
	private static void getExpectedEnergies(StringBuilder buf, boolean doMinimize, String ... aminoAcids) {
		
		SearchProblem search = calcExpectedEmat(doMinimize, aminoAcids);
		
		// print code to make the emat from these energies
		String type = doMinimize ? "Continuous" : "Discrete";
		String resRange = aminoAcids[0] + "to" + aminoAcids[aminoAcids.length - 1];
		buf.append("\n");
		buf.append("\tprivate static EnergyMatrix makeExpectedEmat" + type + resRange + "(SimpleConfSpace confSpace) {\n");
		buf.append("\t\tEnergyMatrix emat = new EnergyMatrix(confSpace);\n");
		buf.append("\t\temat.fill(new double[] {");
		int i = 0;
		for (int pos1=0; pos1<search.confSpace.numPos; pos1++) {
			for (int rc1=0; rc1<search.confSpace.getNumRCsAtPos()[pos1]; rc1++) {
				
				if (i > 0) buf.append(",");
				if (i++ % 8 == 0) buf.append("\n\t\t\t");
				
				buf.append(formatEnergy(search.emat.getOneBody(pos1, rc1)));
				
				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<search.confSpace.getNumRCsAtPos()[pos2]; rc2++) {
						
						if (i > 0) buf.append(",");
						if (i++ % 8 == 0) buf.append("\n\t\t\t");
						
						buf.append(formatEnergy(search.emat.getPairwise(pos1, rc1, pos2, rc2)));
					}
				}
			}
		}
		buf.append("\n");
		buf.append("\t\t});\n");
		buf.append("\t\treturn emat;\n");
		buf.append("\t}\n");
	}
	
	private static void getExpectedEnergy(StringBuilder buf, boolean doMinimize, String aminoAcid) {
		
		buf.append("\n");
		
		for (int i=0; i<4; i++) {
			
			int resNum = 10 + i;
			SearchProblem search = calcExpectedEmat(doMinimize, resNum, aminoAcid);
		
			String type = doMinimize ? "continuous" : "discrete";
			buf.append("\n");
			buf.append("\t@Test\n");
			buf.append("\tpublic void " + type + aminoAcid + "at" + resNum + "() {\n");
			buf.append("\t\tSimpleConfSpace confSpace = makeConfSpace(" + doMinimize + ", " + resNum + ", \"" + aminoAcid + "\");\n");
			buf.append("\t\tassertSingles(confSpace, new double[] {");
				
			assertThat(search.emat.getNumPos(), is(1));
			
			for (int rc1=0; rc1<search.confSpace.getNumRCsAtPos()[0]; rc1++) {
				if (rc1 > 0) buf.append(",");
				if (rc1 % 6 == 0) buf.append("\n\t\t\t");
				buf.append(formatEnergy(search.emat.getOneBody(0, rc1)));
			}
			buf.append("\n\t\t});\n");
			buf.append("\t}");
		}
	}
	
	private static SearchProblem makeSearchProblem(boolean doMinimize, int startResNum, String ... aminoAcids) {
		
		ResidueFlexibility resFlex = new ResidueFlexibility();
		
		int resNum = startResNum;
		for (String aminoAcid : aminoAcids) {
			resFlex.addMutable(Integer.toString(resNum), aminoAcid);
			resNum++;
		}
		
		boolean addWt = false;
		boolean useEpic = false;
		boolean useTupleExpansion = false;
		boolean useEllipses = false;
		boolean useERef = false;
		boolean addResEntropy = false;
		boolean addWtRots = false;
		ArrayList<String[]> moveableStrands = new ArrayList<String[]>();
		ArrayList<String[]> freeBBZones = new ArrayList<String[]>();
		return new SearchProblem(
			"test", "examples/1CC8/1CC8.ss.pdb", 
			resFlex.flexResList, resFlex.allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null,
			false, new ArrayList<>()
		);
	}

	private static SearchProblem calcExpectedEmat(boolean doMinimize, String ... aminoAcids) {
		return calcExpectedEmat(doMinimize, 10, aminoAcids);
	}
	
	private static SearchProblem calcExpectedEmat(boolean doMinimize, int startResNum, String ... aminoAcids) {
		
		// TODO: when the old energy matrix calculator is deprecated, update this tool to use the new one
		
		SearchProblem search = makeSearchProblem(doMinimize, startResNum, aminoAcids);
		
		System.out.println("Calculating expected energies...");
		EnergyMatrixCalculator emcalc = new EnergyMatrixCalculator(search.confSpace, search.shellResidues, search.useERef, search.addResEntropy);
		emcalc.calcPEM();
		search.emat = emcalc.getEMatrix();
		
		return search;
	}

	private static String formatEnergy(double energy) {
		if (energy == Double.POSITIVE_INFINITY) {
			return "Double.POSITIVE_INFINITY";
		} else {
			return String.format("%.14f", energy);
		}
	}
	
	
	// below are auto-generated energies, you probably don't want to edit them by hand

	private static EnergyMatrix makeExpectedEmatDiscreteALAtoVAL(SimpleConfSpace confSpace) {
		EnergyMatrix emat = new EnergyMatrix(confSpace);
		emat.fill(new double[] {
			-11.18489763171806,-12.55046923273086,-1.23183076456355,-13.18241726913394,-1.38102648361967,-12.17157485110978,-1.44054119765207,-13.27642041625077,
			-1.24735585236622,-12.30811570198984,-1.44012131914154,-12.19205283530433,-1.37666346087979,11.97170666214989,-1.00022175683960,11.79027315079589,
			-1.06067460671287,11.64405128556485,-1.01481311321027,12.08349381801739,-1.03286062302938,11.71203358180036,-1.03935939570224,11.60627112765883,
			-0.99997359431684,-14.72861152542798,-1.13416690051677,-15.35780244398886,-0.86790134670789,-14.99101444349932,-1.03913483499934,-15.04432142182633,
			-1.02590361122489,-15.21043488626688,-0.97850203925282,-14.86942266187307,-1.08160753208674,1541248515.54481480000000,-0.24728420809053,-0.29103866137701,
			-0.26126381834802,-0.30141104442893,-0.27931314638851,-0.28231433616307,-0.30255655908424,-0.63193064219834,-0.33030366806118,-0.53891033893614,
			-0.40477988779585,-0.45653629689804,-0.58589428040232,-0.57142702162815,-0.67221439558033,-0.69106927573544,-0.60424457597691,-0.70315957940990,
			-0.63299657849137,141215729.71698140000000,-0.24970844240381,-0.36986883067078,-0.34306845020853,-0.38500053160927,-0.35789708094692,-0.36710420467466,
			-0.38343636968399,-0.67978549686037,-0.42193531294160,-0.57692778144200,-0.48631129345867,-0.51815391699918,-0.61888010538961,-0.66914946551558,
			-0.73492576106762,-0.75865165012662,-0.68427063378405,-0.76187753698787,-0.72047195850119,283.72681273276040,-0.21269621549302,-0.40219057932299,
			-0.35845154903440,-0.39597116388905,-0.38917254524041,-0.37166342220178,-0.40647641184604,-0.69969508946764,-0.45630610194657,-0.60016297990792,
			-0.52123710763255,-0.54595624770364,-0.63684567842985,-0.67179532884595,-0.74644574976694,-0.74538701158534,-0.70236517240788,-0.76027921405136,
			-0.70488472771629,268.86124841467550,-0.21631047168548,-0.41836210320556,-0.38385147224216,-0.41597239707068,-0.40909359851612,-0.39669030333867,
			-0.42288901537605,-0.73214402744585,-0.47074888720532,-0.62687587924474,-0.53990008867930,-0.56574933821199,-0.67014106526352,-0.70222140406418,
			-0.77789945817660,-0.77961746848053,-0.73253873332420,-0.79350275216275,-0.73763634729616,1178912674.16358200000000,-0.25142426010710,-0.29420704276941,
			-0.26013251791457,-0.29567315348006,-0.28329355231338,-0.27603068596580,-0.30118570216769,-0.63574125686637,-0.33446394473995,-0.53668675162898,
			-0.41242972608578,-0.45520225815983,-0.58544421126616,-0.56944541067820,-0.67583704884370,-0.68350198529334,-0.60958346004441,-0.70106119756550,
			-0.62394136048208,112.62535106229947,-0.20407743222087,-0.44970586228839,-0.41333219024799,-0.45507327801232,-0.43667355480759,-0.43365682958134,
			-0.45879279996627,-0.74146935077900,-0.51310703896222,-0.63559810114528,-0.57529253120767,-0.59088457899769,-0.67141050892592,-0.75220278029811,
			-0.80474131211646,-0.81591722968987,-0.76855224520001,-0.82071501263954,-0.78725438778717,114.11927280943505,-0.20073447748328,-0.42688704609526,
			-0.38109405784726,-0.42455003885934,-0.41126580616365,-0.39818839804394,-0.43387348854711,-0.71200080815112,-0.49182138845081,-0.62039932973198,
			-0.54487891503473,-0.57610676519908,-0.64999673623851,-0.70374640091134,-0.76373885344863,-0.76722760315062,-0.72624813087936,-0.77706077155232,
			-0.73492963932901,-28.06028837522226,-0.17187418226261,-0.41696924907316,-0.38606561699200,-0.45113990755049,-0.39206396018205,-0.42900017448562,
			-0.44344020286506,-0.68593162715328,-0.47051428313056,-0.55147580756451,-0.55200242138734,-0.51963492869262,-0.59012763070982,-0.65812318162007,
			-0.67821503007824,-0.71074356032969,-0.65223698189444,-0.70108444124031,-0.69624809153479,-28.70641614639677,-0.13127301528295,-0.43383735410624,
			-0.40109056154799,-0.47406945817337,-0.40480169641659,-0.44970569553919,-0.46478009542118,-0.69856634579333,-0.48686457753157,-0.56559945738824,
			-0.56739055487766,-0.53530145459630,-0.60245260241601,-0.67412761370114,-0.69005462924734,-0.72533373570865,-0.66603583317622,-0.71357651135519,
			-0.71245678819113,-16.70069692439895,-0.15366653794115,-0.51662895152470,-0.44590878602464,-0.53936253130256,-0.47458986503628,-0.49807667342247,
			-0.54263267054113,-0.73566870568414,-0.56316549600427,-0.62219465349343,-0.63109354808779,-0.60519206765037,-0.64214953524984,-0.72588390532742,
			-0.74199571856900,-0.76198597747395,-0.72429265117822,-0.75514270186322,-0.75315373992069,-29.85111670496577,-0.15380238255462,-0.47509504176317,
			-0.43768241266326,-0.50600602893559,-0.44761520402779,-0.48121940154895,-0.50059748061099,-0.73342919016079,-0.52971177220274,-0.59398341526649,
			-0.61265648566245,-0.56944888759722,-0.63011379128155,-0.71224697727130,-0.72569301227435,-0.75174603181630,-0.70740369749032,-0.74226558651530,
			-0.74269587160210,-28.37968845423244,-0.18155700405752,-0.35751074062840,-0.31847723350377,-0.36111138258923,-0.34315370571331,-0.33857826342813,
			-0.36628149450702,-0.62922467857918,-0.40943346403551,-0.50059803493119,-0.48879615973102,-0.46571063498658,-0.53772425827778,-0.57190638789521,
			-0.61451370012747,-0.62296295536646,-0.58513065293244,-0.62692793351249,-0.60042894990164,-30.94635665925974,-0.18485021994442,-0.35307948499364,
			-0.31894397327152,-0.35875420945716,-0.34087429470594,-0.33880156952142,-0.36194676555926,-0.63164184237339,-0.40775315542079,-0.50413181759305,
			-0.48661475668354,-0.46627338187444,-0.54327186366383,-0.57923008230300,-0.62257200878784,-0.63284753766774,-0.59210882391320,-0.63630097648870,
			-0.60934276288832,110.80630737788685,-0.20095350563147,-0.28798412766699,-0.25408402878464,-0.29718354432122,-0.27395418689484,-0.27668739486195,
			-0.29930445444624,-0.57165206117671,-0.32982598995076,-0.44809486709041,-0.41066477956791,-0.40109056899793,-0.48954115555339,-0.49453058090911,
			-0.54948513722750,-0.56384053339501,-0.51006410215981,-0.56835134569886,-0.53294388739645,-30.98033735531527,-0.16976810698912,-0.38721428905072,
			-0.36074200761285,-0.40489905164941,-0.37385673412887,-0.38792352260449,-0.40176124168276,-0.66114709608241,-0.44076210142171,-0.51332287100104,
			-0.52941112295973,-0.48323915163411,-0.55630302772295,-0.62906165991555,-0.65265163229745,-0.67260172732506,-0.63013455678455,-0.66769308646127,
			-0.65928949993725,-31.46316999718023,-0.17342332094016,-0.39091648266962,-0.36604772686386,-0.40335756643361,-0.38126311971809,-0.38828624115618,
			-0.40158456806748,-0.66910631261167,-0.44389235889557,-0.51681288210445,-0.53587868315716,-0.48561574171670,-0.56205927184948,-0.63307463385144,
			-0.65992808287761,-0.67617329865108,-0.63695786441217,-0.67324125506161,-0.66188831493663,-22.12120774427992,-0.17027442935978,-0.38923505029900,
			-0.35371101780010,-0.39679184008795,-0.37552517309101,-0.37585521258174,-0.39936734687380,-0.65525738836212,-0.44533132115049,-0.51836113051749,
			-0.52719123749926,-0.49029835773712,-0.55555028875907,-0.61713093168696,-0.64691552810383,-0.65900984473052,-0.62375623462554,-0.65850495163927,
			-0.64366895670274,8838.47503196292700,-0.20194284106925,-0.29509799245362,-0.26918511324066,-0.31468849503797,-0.28137716543496,-0.29707745249052,
			-0.31109302101711,-0.58213731817975,-0.33901447791704,-0.45324202885173,-0.42150608610470,-0.40721715046233,-0.49697360303982,-0.51982133158361,
			-0.56526065361267,-0.59010249514075,-0.52691109639200,-0.58855568199196,-0.56301576619428,-20.14436389430668,-0.18938659028265,-0.35445502232486,
			-0.33281477579208,-0.37488649482802,-0.34273643379066,-0.35998924646074,-0.36987614310378,-0.63876479605236,-0.40579381824008,-0.49316614420038,
			-0.49527143094977,-0.45618940031881,-0.53908832614183,-0.59856380384536,-0.62908712805439,-0.65331683346472,-0.60009368383453,-0.64871097774130,
			-0.63497208100944,-33.06212155084368,-0.18739132842698,-0.35710009510545,-0.33350854954356,-0.37968663729281,-0.34329446588639,-0.36322449446698,
			-0.37432158863191,-0.63730320991358,-0.40837077749603,-0.49365960990127,-0.49626027745208,-0.45830527062906,-0.53775384288200,-0.60010308912042,
			-0.62796637749433,-0.65354004388087,-0.60018667784445,-0.64773878949980,-0.63650273778071,113.19209269552083,-0.22656404340332,-0.44971543392498,
			-0.42998152415861,-0.47765878409228,-0.43501649124001,-0.46230747501122,-0.46998410321888,-0.73521273315052,-0.49827031520217,-0.59382558180766,
			-0.58414674307333,-0.55461707574879,-0.63955908153160,-0.69895914902960,-0.72886739940383,-0.75949694358382,-0.69726804857110,-0.75252673020207,
			-0.74025424110823,113.70537277011010,-0.23513547302220,-0.46062635272909,-0.44001319316967,-0.48283801069058,-0.44826127196412,-0.46838808630019,
			-0.47706349003590,-0.74898556992802,-0.50813574714633,-0.60537214299186,-0.59623544911613,-0.56477933491348,-0.65261592863983,-0.70835785762439,
			-0.74296923502899,-0.76906259742699,-0.71045803737807,-0.76455578472233,-0.74832009459628,113.04408516686733,-0.23617474396806,-0.42333119794531,
			-0.39456902165862,-0.44163142128344,-0.40780076341556,-0.42274963069119,-0.43933267787735,-0.70590343808341,-0.46984280231885,-0.57964568947831,
			-0.54705711550711,-0.53693082883448,-0.62128339014574,-0.65479630573990,-0.69721286514443,-0.71950098513254,-0.66169279236136,-0.71825635948541,
			-0.69469304188054,95349.16943011780000,-0.26842077800068,-0.33815039046558,-0.32827858483888,-0.37161398659730,-0.32765344549022,-0.36013444571416,
			-0.36000196854770,-0.65289075841840,-0.37352734916064,-0.51822966316474,-0.46446979317867,-0.45781123893273,-0.57121334196918,-0.58307385628283,
			-0.63960849637778,-0.68394368999922,-0.58413861559729,-0.67796876869751,-0.64777622206746,2837.99649210534560,-0.25169548951234,-0.44526858535714,
			-0.43601385642780,-0.47468932657762,-0.43662098198455,-0.46479562687785,-0.46424216341300,-0.74472767954016,-0.49002358392942,-0.59417677594768,
			-0.58405429679162,-0.54875559156470,-0.64671841026046,-0.70064979661299,-0.73492757525316,-0.77017135800724,-0.69812446628824,-0.76276538355642,
			-0.74748596768356,824.69468828884980,-0.23626524373663,-0.43610762655308,-0.42352484937459,-0.46863938886855,-0.42390735058732,-0.45634902851851,
			-0.45779502622953,-0.72951336744870,-0.48150340332385,-0.58317913579748,-0.57235759748495,-0.53985764805963,-0.63289409227490,-0.68975243873003,
			-0.72090479410052,-0.75779742778642,-0.68549502348711,-0.74901808647384,-0.73656715825765,12786046.75555472800000,-0.26328663274568,-0.36579924865450,
			-0.33818778854436,-0.38179729716165,-0.35190728607533,-0.36375551941231,-0.38035119294410,-0.66426203457847,-0.39772213604663,-0.54205156051438,
			-0.48140637242286,-0.48381505639867,-0.58859626968868,-0.57732631487888,-0.64316203133131,-0.66505850148336,-0.59362921185786,-0.66900428764014,
			-0.62688985068062,23849792.74961130300000,-0.27018680824150,-0.35878749134558,-0.33718960713879,-0.37758042821798,-0.34726416712716,-0.36273002021104,
			-0.37366162355984,-0.67054583179080,-0.39129984341936,-0.54677957892127,-0.47485602734173,-0.48399469837749,-0.59699300774134,-0.58903988973084,
			-0.65635297649956,-0.68101254420782,-0.60509737869567,-0.68388674886764,-0.64175768878703,18986109.91319773700000,-0.26379953632937,-0.35012407194349,
			-0.33377681499947,-0.37511757786325,-0.33913506520364,-0.36178268469895,-0.36796944176834,-0.66150671259894,-0.38776651926993,-0.52713110889255,
			-0.47930047017550,-0.46818762666919,-0.57896372743883,-0.57915587163338,-0.63888661148079,-0.67203341110888,-0.58719101895297,-0.67086707112839,
			-0.63543325495734,5682.62148228075200,-0.21837001972242,-0.53115797047149,-0.51689425371324,-0.57877024916031,-0.50948295572377,-0.56247459199639,
			-0.56381746108676,-0.81805732240286,-0.58135934809279,-0.68021544025128,-0.66576157254043,-0.63960143065517,-0.72546821310134,-0.78880562353933,
			-0.80997765103327,-0.85736250272300,-0.77612886595216,-0.84348022251510,-0.83860767759671,1726.80340673180880,-0.18966271538106,-0.54649252277207,
			-0.53002982703301,-0.59947360040555,-0.52055618744278,-0.58110699841807,-0.58316969037874,-0.82844509423349,-0.59625460273872,-0.69155769768753,
			-0.67979388468272,-0.65299529827772,-0.73505747375287,-0.80268654975524,-0.81964230935799,-0.86906537351389,-0.78799241765540,-0.85326280753329,
			-0.85226341656112,2580.55801061422430,-0.19331316981355,-0.49633080066153,-0.48656447325351,-0.59622967120688,-0.44758893830255,-0.57003584859166,
			-0.56621526143509,-0.77615376264923,-0.54822732011801,-0.67227463143059,-0.61248540596434,-0.62622801888351,-0.70743621937405,-0.76055774559167,
			-0.76635340685074,-0.85682659524059,-0.72174257804727,-0.82807751738574,-0.83416324317002,-29.25215096964729,-0.19746729933866,-0.56978477659489,
			-0.54156022492344,-0.61495293531069,-0.53964034854220,-0.59209592173620,-0.60367722266060,-0.84015082399944,-0.62193500565843,-0.70904717598058,
			-0.69989529377376,-0.67571367285160,-0.74752939115246,-0.81798684366089,-0.83543435591862,-0.87424398665082,-0.80863483101540,-0.86177514916392,
			-0.85972523248041,-29.39119781728211,-0.23589716293232,-0.55259024594659,-0.52528390393032,-0.58970772893807,-0.52786512693435,-0.56910724791281,
			-0.58043534705009,-0.82863163146993,-0.60330094750724,-0.69568501650672,-0.68350692580148,-0.65931358437132,-0.73657381864880,-0.80063607998838,
			-0.82494680538788,-0.86102131875635,-0.79437173514244,-0.85107071862772,-0.84362953307968,81.89245964750002,-0.03182133111382,-0.74747421617025,
			-0.75315273256034,-0.74364004869433,-0.74632344222724,-0.74502379752850,-0.74695884455207,-1.05627491933054,-1.02511662237201,-1.01136329740440,
			-1.04602478295924,-1.01579404537161,-1.02042613641896,-0.91953488791772,-0.91755137743839,-0.91542918478265,-0.91921177046283,-0.91594192633398,
			-0.91734887978269,-1.54672409227993,-1.54818317732530,-1.47879405767539,-1.54322105459145,-1.53688566730276,-1.51907739732100,-1.45125752551442,
			-1.74500801533955,-1.72523914525262,-1.75493563960914,-1.70568269310154,-1.55557838577465,-1.64436532501969,-1.63798612158748,-1.61103745878478,
			-1.62643861526652,-1.56194940399688,-1.65956769308923,-1.65610434509466,-1.63066817308254,-1.43503090478411,-1.44654263184984,-1.41256118751860,
			-1.46108754372407,-1.45515963042709,-1.44853745590101,-1.43119342195865,-1.46186909937779,-1.46163760792707,-1.47492997381811,-1.45153909794228,
			-1.45050177422022,-1.42565094515494,-1.46068834592261,42.67793611075126,-0.03670613008751,-0.64283427705032,-0.58918414548732,-0.64187350450898,
			-0.62751039474925,-0.63343464228955,-0.64490901860713,2.64460320929424,2.65784434000583,2.66894909800905,2.64238827973618,2.66548667326362,
			2.66450923133723,-0.76839215103904,-0.76617491216238,-0.76473300384984,-0.76786799655385,-0.76477522119119,-0.76657827440660,-1.73671356139369,
			-1.74332563448015,-1.67227887507502,-1.73925329845380,-1.72754519158589,-1.71442649404688,-1.63673288238469,-1.91397124824572,-1.89872533431988,
			-1.89772455667871,-1.87730527860897,-1.72332067327840,-1.82061868375224,-1.80258593675884,-1.79138248971414,-1.80079037009493,-1.73590111721330,
			-1.83527780376673,-1.83109748202434,-1.81127466194518,-1.61144838588215,-1.62388430108068,-1.58293209421092,-1.63909239902928,-1.63284884955531,
			-1.62559735146243,-1.60281695750682,-1.64004682147612,-1.63983951474981,-1.64918527055814,-1.63433824402839,-1.62565696037493,-1.61241486459721,
			-1.63428630357045,38.48377127974877,-0.04524725550607,0.05375117984924,0.09227195782867,0.05200046375443,0.06754832322461,0.04459851445438,
			0.05159022539897,4.84774201500523,4.86932522955307,4.88061476347180,4.85184643015874,4.87747688205371,4.87139487218067,-0.17365265915241,
			-0.17156018802433,-0.17012273030389,-0.17314908900439,-0.17028348350786,-0.17183151678597,-1.89419816796433,-1.90068328466491,-1.82888136213901,
			-1.89552890448112,-1.88403616191302,-1.87119172138697,-1.79957763014750,-2.09252770108038,-2.06482695540870,-2.10155105132564,-2.05013402654071,
			-1.88860899481071,-1.98633248794519,-1.96859918765707,-1.95204470906066,-1.96716670711963,-1.90177951736636,-2.00240069074560,-1.99764701769847,
			-1.97246124195639,-1.77286708187271,-1.78405117115679,-1.74330792369627,-1.79968406086091,-1.79297512632644,-1.78687827515947,-1.76302915324872,
			-1.80005385123951,-1.80034301892411,-1.81711179070784,-1.79425886288441,-1.78336248822865,-1.75997284355905,-1.80267602444759
		});
		return emat;
	}

	private static EnergyMatrix makeExpectedEmatDiscreteGLYtoGLU(SimpleConfSpace confSpace) {
		EnergyMatrix emat = new EnergyMatrix(confSpace);
		emat.fill(new double[] {
			-10.65318811007568,-9.48840986632377,-0.96373042704017,-10.19699962234164,-1.10792120662228,-9.06355106032938,-1.16986344767586,-10.29478723244586,
			-0.97453888706844,-9.24630018155107,-1.16804149726335,-9.04352335562025,-1.10838036932644,-12.83574477814928,-0.95523306990956,-13.08023833775742,
			-1.01771752577905,-13.13496232521485,-0.97233574373401,-12.75439773179176,-0.98809961056195,-13.10938895040057,-0.99749157373287,-13.14377785915285,
			-0.95617923755029,-12.91361782276331,-1.21408766586040,-13.59693974330352,-0.95071024997953,-13.12671775034629,-1.11562457238726,-13.29547289661611,
			-1.10849274080255,-13.39513386701661,-1.05805032685466,-12.96767994724094,-1.15816321122515,5.07372688576625,-0.15071016134613,-1.55530846501363,
			-1.65987654199462,-1.58315300212243,-1.61271502604976,-1.64462323616551,-1.54712736441454,-1.87604009071651,-1.49794738078100,-1.30008603988581,
			-1.93225793044781,-1.31025974867988,-1.46050103724922,-1.64629781125589,-1.58549268202013,-1.61636648513700,-1.60772139986169,-1.58700223535503,
			-1.65434864480859,2441.06984573589300,-0.16790301941165,-1.35798864066110,-1.45703623348021,-1.38905517190382,-1.41009115347261,-1.44718724160105,
			-1.35319412226969,-1.68517512546167,-1.26110023091889,-1.19946604818975,-1.59733539109005,-1.17317350638913,-1.35511671116641,-1.43440761233292,
			-1.40742468024897,-1.43977557616905,-1.41060041127574,-1.41634516681838,-1.45742511896669,-19.00367517068371,-0.16901097185870,-1.15166751586884,
			-1.22171638222961,-1.16108934620205,-1.19844254428638,-1.20145822060043,-1.14068758934432,-1.42203497253102,-1.07934225913293,-1.10699917554447,
			-1.28889391315092,-1.06889922697633,-1.21562224562879,-1.19595387907828,-1.20180870161851,-1.21691495483254,-1.19312012956669,-1.20828469306057,
			-1.21617189080760,-19.15946166389191,-0.16983607546848,-1.15166271621462,-1.21019633874282,-1.15266511513130,-1.19677229519761,-1.18735215210646,
			-1.13747161253813,-1.40958091932939,-1.08245718338660,-1.11152363428979,-1.28196645120803,-1.07567326790489,-1.21201839038301,-1.18689939002075,
			-1.19723555252369,-1.20638531873142,-1.18829361040608,-1.20098053841901,-1.20423516361887,73.96971843893068,-0.18571023760108,-1.34796484546360,
			-1.39702015015070,-1.33003104301890,-1.39949853583674,-1.36356801925040,-1.32170487192454,-1.59585887440529,-1.27663657884945,-1.28801279287249,
			-1.48289577480364,-1.25846802265174,-1.38764619227744,-1.37777259409524,-1.38865475597069,-1.38196916078091,-1.38747416810599,-1.38107240479322,
			-1.38364792762667,67.12739922315822,-0.20921577113230,-1.19541077095022,-1.23598230780286,-1.19419643695493,-1.23168231026597,-1.21987153565035,
			-1.18274913064038,-1.44680828181270,-1.11438384714004,-1.17813510264184,-1.29752255065425,-1.13398947364130,-1.27071835878509,-1.21323193315081,
			-1.24104234700560,-1.24991918332585,-1.22161932754118,-1.24761170568502,-1.23903119637744,3.96252445467034,-0.11247475217984,-1.35692172655284,
			-1.38968679151430,-1.36844394551659,-1.38028106283962,-1.38906286083755,-1.35401831312890,-1.56647340184855,-1.29554692368808,-1.27667392266536,
			-1.48595812649238,-1.26267393516527,-1.36211223841447,-1.38603206852716,-1.36775965685152,-1.37640438044037,-1.37794616150592,-1.36497459191727,
			-1.38966453505096,5245.79060070078300,-0.10182981018558,-0.15136148162853,-0.40709764664771,-0.16619058311038,-0.29915540678933,-0.32559504339188,
			-0.09603187574702,-0.46014399980441,-0.21739968309981,-0.30933363097626,-0.31892819265588,-0.21751240097495,-0.44196725649488,-0.33384038152299,
			-0.38543533656049,-0.39604048816590,-0.34640626917434,-0.41493960066801,-0.35451732279611,-1.29390603351609,-1.28697487359207,-1.06368443848994,
			-1.43246794091368,-1.31061868487753,-1.35334460233142,-1.44092711500881,1726.26099454015000,-0.18605123206791,375.15251509897900,374.02874454823944,
			375.40582695968925,374.36166369315270,374.93687321989535,375.45479395002740,374.46026963439573,375.00936135206194,374.97668913044834,374.74328037842950,
			375.07971833606615,374.68066076517210,374.72646625529455,374.66684241649676,374.71223148254006,374.69448503100200,374.67650001674100,374.73480255481815,
			-1.97937329373590,-1.99049313754835,-1.67645979795281,-2.01987761986692,-1.89824988989542,-2.04382468808382,-2.33964136751240,2718.05559694319600,
			-0.04950387667063,-0.23266927082478,-0.56810753612186,-0.40470195326073,-0.34218008960436,-0.61839609573843,-0.25323855334192,-0.58850746178095,
			-0.15262190486306,-0.33431735754077,-0.31258738816708,-0.16195814657902,-0.60280272377129,-0.47347778199102,-0.44387904896387,-0.52066055290298,
			-0.44151904821821,-0.47975958957544,-0.51637719953221,-0.63021053969031,-0.54110148103653,-0.50855053714151,-0.96805044015998,-0.87279924430988,
			-0.90248091587926,-0.90995096801178,1902.21562131460810,-0.07437980577009,-0.05524837761029,-0.38297035416083,-0.19440814517831,-0.18139299601489,
			-0.37999972403214,-0.06241131478745,-0.40458861165204,-0.15121840315425,-0.27247935527457,-0.24882676198704,-0.16279985571396,-0.41222613283014,
			-0.29656174957984,-0.32737043942206,-0.36831815351672,-0.29569495893651,-0.36854805243378,-0.32805191000747,-1.12478916986825,-1.10881690765437,
			-0.91343338581075,-1.30301022285698,-1.17924997515422,-1.21090996531564,-1.25419746021823,545802.31291239660000,-0.05909783446426,-0.27893611652216,
			-0.46112924697147,-0.30162044074515,-0.36671329978481,-0.40383943310574,-0.25694008381973,-0.46697812629019,-0.28108751632866,-0.31291468954978,
			-0.37057328105093,-0.26493959533331,-0.41011802878509,-0.38506796391891,-0.39055913705311,-0.40264511092315,-0.38304671690233,-0.40108745743253,
			-0.39407299786222,-1.21838406986246,-1.18260723727042,-0.87211887531145,-1.28168937343322,-1.23983806341063,-1.25077994967529,-1.34081769293345,
			-15.45205898009560,-0.25859324316813,12056.72483063713800,12003.17664924472400,12054.86644715198000,12052.55718474262200,12033.69626262292300,12057.15385396645800,
			2151.89869589134470,2154.19154010692400,2153.32996237760740,2153.55597344700440,2154.15128257793460,2151.36488263923770,1897.32020923900700,1897.17438568158560,
			1897.05072952598830,1897.30172485828370,1897.03275418279700,1897.21168204596100,-1.75916725464850,-1.82370845127459,-1.65278980893700,-2.02842643350704,
			-1.85155816652490,-2.03112167252318,-2.17056363671770,-16.54627589422750,-0.08762901121350,0.15983233066281,-0.44811903215882,-1.03160653262812,
			0.05025395637943,-1.38465493895178,-0.13909785956468,-0.71177293914129,-0.23841813767934,-0.73286587282225,-0.36417284830438,-0.34893739029245,
			-1.02343094831565,-0.49225070047615,-0.56193241652162,-0.73972346569945,-0.47893443222361,-0.70763925245169,-0.60095103300691,-1.30962268216402,
			-1.31294196711988,-1.15665994281304,-1.54260534726097,-1.39944737794549,-1.45725725930083,-1.48757875674023,127.25515680955156,-0.05787641189962,
			-0.38567162673612,-0.98720655630115,-0.76781899275073,-0.55147976484023,-1.18259649722455,-0.44996907705771,-1.13139476307692,0.10106267744621,
			-0.62204591631249,-0.23133451801997,-0.00201629640617,-1.63371472016668,-0.74013263677846,-0.65243933128166,-0.83673200950136,-0.65741067264030,
			-0.72584344453317,-0.84776294286861,-0.60773998458124,-0.55651983923336,-0.70256524991476,-1.15821306579854,-1.03558566790476,-1.09592903254503,
			-1.06429867797952
		});
		return emat;
	}

	private static EnergyMatrix makeExpectedEmatDiscreteASPtoPHE(SimpleConfSpace confSpace) {
		EnergyMatrix emat = new EnergyMatrix(confSpace);
		emat.fill(new double[] {
			136.91563957788102,559.29790051081190,-20.01154802696374,-20.01198042216748,-2.39669908910579,42.82594134270074,-1.21295088362082,-1.23084373713685,
			-0.95560046360195,-1.20573622005578,-1.17058571632116,-12.39150992112944,-1.05125209762279,-1.06658788521491,-0.89463400458684,-1.11414327070703,
			-1.07701309275459,-12.79143450907142,-1.06373082975283,-1.08307982405098,-0.99779196771016,-1.22505448233685,-1.15460441738074,408224729.43032587000000,
			-0.23900191769386,-0.24479061898893,-0.20185518858011,-0.22634747557822,-0.24031577730424,-1.28887727471182,-1.41390509731185,-1.28278045740886,
			52.43275582969316,-0.16226778131892,-0.15854005064979,-0.13910023552726,-0.15613989704206,-0.18690436154100,-1.29162882005074,-1.40791832122250,
			-1.27782108160934,-10.92218769395612,-0.20746632613133,-0.20144009071981,-0.10636117701114,-0.14099940154363,-0.20547610826489,-1.18795202990427,
			-1.26025537379859,-1.18392396576569,-10.07337523535645,-0.15396071781901,-0.14852160117731,-0.03304204643123,-0.06120931521731,-0.16513518413150,
			-1.18629181245787,-1.25672529385779,-1.18049532205225,-11.59136935985416,-0.17612160942608,-0.17526294322529,-0.13460886703984,-0.15286603597355,
			-0.18620359759148,-1.17579986214757,-1.25033998708713,-1.16301840173653,-11.03776243931168,-0.16035310660030,-0.15796146648929,-0.12866255075544,
			-0.14678417624791,-0.17865471530706,-1.18264496780689,-1.24479537705922,-1.16866965421061,-11.32320022861451,-0.21120638866160,-0.20353679262588,
			-0.17784888517037,-0.20113644825798,-0.21076768711479,-1.18178889830514,-1.24576686701627,-1.16921499292759,383.23866566955790,-0.17106373652932,
			-0.17675820485250,-0.16149025261347,-0.17561279903910,-0.21819740990022,-1.29111650267352,-1.35584136474421,-1.27710649573492,800.98544671315810,
			-0.24267762071892,-0.24752089164825,-0.24108534894282,-0.27665048356172,-0.27076333721402,-1.28850510564971,-1.35335287154758,-1.27848094627031,
			67919.35227738268000,-0.22773231439023,-0.23866604709810,-0.19668739460624,-0.21579084705481,-0.25685506351328,-1.28273692558160,-1.36494099419029,
			-1.27636739397819,-11.90791063302230,-0.41605363252468,-0.42574665306794,-0.88845896040935,-0.90083390430346,-0.41933268015104,-1.64457470774588,
			-1.70994441834337,-1.62820905878600,-9.49215552365857,-0.27873095981650,-0.27847752402723,-0.20843857419156,-0.31099522747073,-0.27829988721827,
			-1.34583469152906,-1.42073397972546,-1.34515292186495,-11.31469624039974,-0.16801355501925,-0.16636728737280,0.01976582911763,-0.00525256720574,
			-0.20020853096790,-1.34069051577437,-1.41699220856336,-1.33626030499713,1288151813489.19950000000000,-0.08636171866473,-0.08367056834729,-0.07244318989646,
			-0.07099167586814,-0.08211501909635,8.16731436170320,8.45218224132250,8.54847647197311,-1.48112813298378,-1.49788852574916,-1.77491981281187,
			-1.84258401495949,-1.62391518134033,-1.62195176214952,-1.61250520153696,-1.41691718642638,-1.40568194961534,-1.40946869009475,-1.53991728395438,
			-1.47340469352516,-1.50470752344386,793844508.97293730000000,-0.03524242498045,-0.03542960285622,-0.03289561512077,-0.03696833390711,-0.03199417986834,
			-0.76753086172962,-0.70000910029204,-0.60542140961335,-1.19973869211338,-1.21663869349803,-1.45816516439999,-1.51488016769932,-1.32692415243220,
			-1.32113364242567,-1.31124313494756,-1.12127476676849,-1.11325344223567,-1.11809013508423,-1.18311263511513,-1.16522268356357,-1.18679809020546,
			83.94087533735910,-0.12813935538281,-0.12579199588672,-0.11835196686199,-0.12580567495505,-0.13363660086842,1622075917587500.50000000000000,49006.27151701816000,
			61566.18451647393500,-1.76533710689798,-1.78131652809427,-2.00000487185078,-2.05241361919381,-1.86577209033048,-1.85239400517333,-1.84546948304036,
			-1.65668651527701,-1.64780670240117,-1.65852219857688,-1.74905589898364,-1.70975955443977,-1.73233917124751,56.15946609712546,-0.09499178450693,
			-0.08858475504412,-0.04645827315787,-0.06906335114313,-0.10146667270866,761.97913898562830,61.41020703329888,-1.32164233370521,-1.68439452912058,
			-1.69880305297047,-1.89325498422944,-1.94482661129807,-1.76657039980421,-1.75497118991405,-1.74906249951254,-1.55383662209031,-1.54601354897387,
			-1.55487680097992,-1.63552143985523,-1.59707482627462,-1.61836980579127
		});
		return emat;
	}

	private static EnergyMatrix makeExpectedEmatDiscreteGLNtoTRP(SimpleConfSpace confSpace) {
		EnergyMatrix emat = new EnergyMatrix(confSpace);
		emat.fill(new double[] {
			679.00187379223470,237058324.33089966000000,243.60747236742307,69.18531107973202,-22.48049509545349,157.34498616854620,-22.78029633171821,61271.30319413747000,
			21289.07875802994400,953.73425237212780,-1.67445159407726,-1.71573890650799,-1.94518893730012,-2.10529658390733,-1.79398313744014,-1.57881244534732,
			-1.58064530214799,-1.60000708821825,-1.63043601966865,42313.16372166756600,-1.67166480507315,-1.71292717609072,-1.93525240422886,-2.09408934457454,
			-1.78727945351755,-1.57677852602655,-1.57818018842592,-1.59765834566852,-1.62749299834605,44.18538688010391,-1.38690469825881,-1.41909032311011,
			-1.64859962704689,-1.81253061549784,-1.51003824962536,-1.29204043849390,-1.29806354625158,-1.31286011841921,-1.34336350189534,83.72908372336694,
			-1.32127058742286,-1.35179968046621,-1.58623687015264,-1.74938358907646,-1.44707332929337,-1.23332854935390,-1.23928791964156,-1.25405008688353,
			-1.28408176073245,466.66825498616130,-1.89123470236094,-1.94067302390542,-2.12809424238253,-2.29067535326726,-1.99171936781707,-1.79138581205503,
			-1.78850898658003,-1.80807506133657,-1.83773466174121,197.74910497190209,-1.52473302550821,-1.56778322726136,-1.76153999230709,-1.92307386474604,
			-1.63009335527106,-1.42634027234681,-1.42692222636921,-1.44087792458959,-1.47044372833474,212.73080208987963,-1.56653171377788,-1.62356131659509,
			-1.77953670339314,-1.94140856875889,-1.64748247319770,-1.43691993178710,-1.43689236970515,-1.45128414936636,-1.48160833797491,3428848.66410222040000,
			-0.22337552889401,-0.21951701785512,-0.25967892285254,-0.27968361949194,-0.24981334830609,-0.20568613006600,-0.21225532479666,-0.21629890233305,
			-0.20675313231760,-1.82027912970149,-1.53942581294993,-2.22423863228186,-1.98678470764122,-1.80948564651063,-1.81230234017135,-1.83633411094727,
			-8.66089448291007,-0.21376294992652,-0.20782428678460,-0.26004301290679,-0.27865647814314,-0.24734377689664,-0.19559714271690,-0.20328105542808,
			-0.20887136985977,-0.19903486531051,-1.63760734856773,-1.37632689183755,-1.79041783376525,-1.69526659826931,-1.53145449170162,-1.53180252792169,
			-1.54253554021512,1.69393071957849,-0.20888013321921,-0.20294480004327,-0.24968964311577,-0.26751290639034,-0.23632825373951,-0.19114513908744,
			-0.19874540312768,-0.20428105968349,-0.19328234729667,-1.61272247616934,-1.35098859813501,-1.77080214226021,-1.67431989505106,-1.51068251773736,
			-1.51123003986093,-1.52144868803929,32.99802348550968,-0.38388579521989,-0.37749135008168,-0.48879036385502,-0.49450292256978,-0.47931003674338,
			-0.32633465833042,-0.33766480818397,-0.34407838844837,-0.33058630459205,-1.95476634910136,-1.69330812768007,-2.15290403823797,-2.02291990593020,
			-1.86349965340800,-1.86123192234758,-1.89067831156045,383.10659717334320,-0.29471315372792,-0.28702780108231,-0.38180246612821,-0.38993253088541,
			-0.36814916079954,-0.25885677059906,-0.27012046539364,-0.27553360929039,-0.26274361093852,-1.80873042417236,-1.54879373064681,-1.98973830197891,
			-1.87377281132181,-1.71200958409406,-1.71088597344658,-1.73191760046334,717036899.47906780000000,-0.04516327952719,-0.04487134676641,-0.19923884254955,
			-0.02143197016220,-0.04716392010192,-0.03982952339505,-0.03904322568000,-0.09401397355316,-0.06494097961007,1.54162476889448,0.35479760551517,
			0.37968748565918,0.34802322358581,0.69355577193747,0.69371914043255,0.74992968000720,-1.68393844202200,-1.73304395730621,-1.79734498363826,
			-1.41833635742175,-1.43509420429535,22725408723.76246600000000,-0.32614467876031,-0.32809760027675,-0.60922411217308,-0.04730370772798,-0.46926333945308,
			-0.33794297725717,-0.34137151975930,-0.48794231091673,-0.42717179404416,3207.36837126754650,3165.73379943441020,3156.22640703882050,3156.18648578284000,
			3135.26951137064400,3135.31017587087760,3135.39682554156750,-2.17103156809522,-2.24774909222277,-2.28734887559772,-1.95191590257144,-1.97559134331393,
			1375700119120607490.00000000000000,-0.04035023172643,-0.03995777218990,-0.06535448104951,-0.05334100000439,-0.04630081848956,-0.03833555198649,-0.03943813482959,
			-0.04612666080842,-0.04340330812946,0.15140459322352,-1.13127743235901,-1.11096107244477,-1.19697286958191,-0.91418137971545,-0.92690033330756,
			-0.85742867203200,-1.56185791913413,-1.53671786312096,-1.63020685197847,-1.24708368095457,-1.25767665974279,438974350.63024080000000,-0.03827489160543,
			-0.03900543078361,-0.06398359744351,-0.05435477420704,-0.04702284634958,-0.03924145607336,-0.03974689092332,-0.04429592837589,-0.04346857084653,
			0.33666307604013,-1.01295331478709,-0.99149451319670,-1.03484508328287,-0.81262597359562,-0.81977089263553,-0.76486801599999,-1.45710286611364,
			-1.46190694753090,-1.54933833602119,-1.17130337658125,-1.18227818779303,114.62396641251748,-0.76180017464996,-0.77939351319932,-1.10117689969265,
			-0.97154725571787,-0.83252274764881,-0.76274278836285,-0.75474799551149,-0.80769943860626,-0.79739210745785,748205963.64450840000000,147319234.84068840000000,
			736564596.99923260000000,736564678.34192060000000,9432787.80529877500000,9432840.50299201500000,1973090.18087962850000,-2.49264447201335,-2.43412927976500,-2.51067687826100,
			-2.19197437012941,-2.18638131543865,201568.01455711323000,-0.12828296220590,-0.13338925498460,-0.15195258154029,-0.17592150127811,-0.14948369680700,
			-0.13751703553499,-0.13552473900744,-0.14035891545857,-0.13969730643010,1342632156.76014760000000,330798.68925840340000,172.49579216001231,206.74765978438400,
			267.56521423411570,267.50300446883375,50.26094102192376,1.45562612056979,1.55179218418795,1.45707208055091,1.83596885295037,1.82835091853957,
			20476946145.63862000000000,-0.05630288525840,-0.06663929172774,-0.14503517506108,-0.05987836049370,-0.06308537272916,-0.06103704169446,-0.05422909317022,
			-0.07978393299106,-0.07799353809910,112527320.48129016000000,1053869.02071038600000,38216.11735476660000,38304.21119773361000,41.83547994879611,41.52924814936535,
			36.16583026158611,-2.04804024950477,-1.97577229090389,-2.06812397718917,-1.68924828368762,-1.69820647267188
		});
		return emat;
	}

	private static EnergyMatrix makeExpectedEmatDiscreteHIStoTYR(SimpleConfSpace confSpace) {
		EnergyMatrix emat = new EnergyMatrix(confSpace);
		emat.fill(new double[] {
			16721948.50390606200000,9110858.17040926400000,-11.56113038428850,-12.46904571981044,-3.03903472034233,132.94497022015870,93.33673779249861,-15.27124703159180,
			2144967.85985721800000,-0.05774946397647,-0.05776996797113,-0.32279357628002,-0.20872102244736,-0.17086536258680,-0.02794452398232,-0.04532261270417,
			-0.04854477688150,129183009.37787631000000,-0.04383980444833,-0.02064874717327,-0.14097204476651,-0.07106992945482,0.00557101965198,0.05341010795957,
			0.00927430303068,-0.00487113887105,70910.54210738631000,-0.14659285385483,-0.13819217374773,-0.34501945986342,-0.24502475400303,-0.19045413807180,
			-0.07124642727720,-0.09890698369692,-0.11119099083336,45061.55471765097000,0.14042884936833,0.14738441401117,0.02767974507335,0.09314873194500,
			0.20208419001681,0.29778917702476,0.26247381761052,0.22775190583203,53162.63523389290000,-0.17207893626660,-0.15295433330520,-0.39632333839401,
			-0.27863699471334,-0.25532189327348,-0.14785133422600,-0.16842691076673,-0.17414837300121,49225.41363837266500,-0.05995383216587,0.03857455121971,
			-0.30483025686143,-0.10587445503919,-0.19561033219511,-0.16739115023030,-0.17607946025246,-0.18098541625631,11463.81887068581500,-0.07247519691977,
			0.01805558448696,-0.41348862410104,-0.22076384969822,-0.29293655363231,-0.22798041064903,-0.23050486958401,-0.23120865294975,3062.89884331063830,
			-0.12424203423978,-0.06437917825023,-0.40377816341702,-0.24155495248196,-0.28280434365888,-0.20088655734831,-0.20786307194703,-0.20902152464461,
			8986.22243023858800,0.00001262787255,0.03865378867819,-0.15023447828485,-0.04476329835221,-0.03719572040984,0.01972027809545,-0.00403416658434,
			-0.01162412216287,606854.24071167970000,-0.06664773241584,-0.04601183973554,-0.27931410688169,-0.16777777102850,-0.15090709572771,-0.04802947528108,
			-0.06464718177274,-0.07060887311753,3825.52786938312700,-0.05806695530361,-0.01243983375552,-0.32982317563467,-0.18274642040167,-0.20822074244817,
			-0.12210674859380,-0.13010373793968,-0.13184003336011,171128351.76079917000000,-0.14366579488061,-0.11877249342839,-0.39006155151460,-0.26550638980138,
			-0.26604210762913,-0.16065911953208,-0.17028163903581,-0.17379666012289,1783.55243869672000,-0.08444148503272,-0.04277051112298,-0.29305558982842,
			-0.15952277465582,-0.17245992982952,-0.10338956080904,-0.11716178728179,-0.12273408090371,493.14913810301430,-0.03711322011003,-0.02040068737417,
			-0.32108004084979,-0.19482936356652,-0.19414295818640,-0.07070515453242,-0.07931262299487,-0.07673834916723,779.17273257808070,0.08541073307345,
			0.12787967355051,-0.25488800937660,-0.10175336249285,-0.13146870441314,-0.02832452902068,-0.03520588652665,-0.02835742761199,4072121.91691076000000,
			3.31899292940111,3.13800169825516,2.48533934285352,2.66119884432362,2.61578133755707,2.85835217938846,2.85618977217945,2.99848394845475,
			1539.84258600346540,-0.23609384645623,-0.19677185804939,-0.58090933210592,-0.42260070892813,-0.45422475213771,-0.34725876587858,-0.35536611477861,
			-0.34493242983211,898.18237129514720,-0.37302355550952,-0.35746395368458,-0.65280467870083,-0.52343169292918,-0.52284539613046,-0.39852433068674,
			-0.40862361793461,-0.40383590151621,1659.79587562189770,0.08644276793374,0.03088012269133,-0.37440357011050,-0.25042611155863,-0.24186837711810,
			-0.04579023079555,-0.04517599930307,0.00462608837723,500373281.15036360000000,-0.30624945169420,-0.31939871723719,-0.64591308694281,-0.52004183798258,
			-0.51217352161673,-0.35426780497286,-0.36009721792215,-0.34248328505711,20097.45763274872700,-0.11264504602166,-0.11424990207753,-0.56166095465936,
			-0.40939737618155,-0.43222289760966,-0.29223010458792,-0.29671676451150,-0.26769335273312,81838.74690774774000,-0.27520667283050,-0.27822278829432,
			-0.53806475005633,-0.42561594875546,-0.39562710999312,-0.25161095837910,-0.26516710884596,-0.26366895621785,629714.52852543700000,0.13087396542040,
			0.10252268368973,-0.17821912536863,-0.07479706614443,-0.01713958382781,0.15573096461189,0.14744939556527,0.15100329766620,12803.76510144592800,
			-0.37015136652115,-0.36767846778388,-0.84393303447035,-0.67470328597847,-0.71482790969725,-0.56840004649900,-0.57238281707514,-0.53451479424266,
			1032827.35835145170000,-0.27334920163368,-0.23987060766009,-0.89956784050497,-0.70591251334247,-0.77501778159964,-0.65649572086467,-0.65629411580436,
			-0.61002096271787,23829.78869004675700,0.09158785883334,0.01116548925781,-0.71982486752254,-0.51306589961533,-0.59190615157105,-0.41748431348854,
			-0.41801522171570,-0.30966525813851,9346.95456429938200,-0.44894942578867,-0.36340741005021,-0.87620543164821,-0.67146146802752,-0.75149020743025,
			-0.66845055020353,-0.67209694697974,-0.65824193435528,36641398.54575956600000,-0.25383231069734,-0.23263540925137,-0.27956815783606,-0.27186075681720,
			-0.26031254032015,-0.23477879021444,-0.23486874049462,-0.23593237286196,-2.38933101852878,-2.25589037090300,-2.36940813401995,-2.34556605930983,
			-2.43667960804525,-1.68542751916117,-4.12722646107017,-4.10257364545811,-2.99866861179575,-2.88895543041819,-3.03387520379428,-3.01395329205776,
			-3.52688763543883,-2.86024105114492,-2.79432133778847,-2.28487846517793,-2.24690760859795,-2.31476753079852,-2.27633085247480,-2.29968541219706,
			-2.23091628120669,-2.28515593424478,-2.26963012299484,-2.32416871083522,-2.24052002073909,-2.28303018567997,-2.31381117776222,36641400.26061017000000,
			-0.22787235055266,-0.19958322572234,-0.27038969350975,-0.25564717592730,-0.25100227711750,-0.22898562505818,-0.22810015832274,-0.22887339194832,
			-2.35951449093924,-2.24274333513428,-2.35141645324004,-2.32348469480219,-2.42089881300294,-1.53696301309028,-3.27422816663068,-3.92222947143223,
			-2.99965555751289,-2.87149625081441,-2.90704773404202,-2.99604489501329,-3.56432706878544,-2.80904215657365,-2.66931686618820,-2.22025313202415,
			-2.12485750145928,-2.26264761634408,-2.22451717554892,-2.24952463820449,-2.14598134797526,-2.25045209552275,-2.23044921732356,-2.23554225692279,
			-2.10273348157130,-2.19770874299573,-2.02326277724063,394.39194433940054,-0.19358824718919,-0.17953531473212,-0.20218263676721,-0.20227557148979,
			-0.18515802751584,-0.18112515151646,-0.18205602295803,-0.18300698436704,-2.00031696094399,-1.90383888178846,-1.99242745479679,-1.96792798555084,
			-2.08429214546553,-3.22938422797145,-2.39485607284764,-2.52628714164949,-2.43438480522675,-2.30459369456547,-2.38537213314823,-2.38998814746292,
			-2.90884755545142,-2.25769840013255,-2.23076227029782,-1.86659023115032,-1.89140688536438,-1.92326293918189,-1.86287634027211,-1.88352465458477,
			-1.85700612084459,-1.87513244538235,-1.86027504797209,-1.87922705347263,-1.84331499061950,-1.84424381507057,-1.86557932702620,394.68215544177264,
			-0.20496605505131,-0.19210249626062,-0.23705649249085,-0.23796511003286,-0.21984160074001,-0.19429486812796,-0.19430539648297,-0.19498556682310,
			-1.99448999039182,-1.89742157170738,-1.98607828819006,-1.96544406155165,-2.07271852516991,-3.19949743685539,-2.37685602272904,-2.50751676017560,
			-2.41561448359200,-2.29205474331814,-2.36664220142335,-2.37279760293700,-2.88164376748981,-2.24757313954600,-2.21982046043772,-1.86443699306660,
			-1.88140822913735,-1.91347950923871,-1.85906480848711,-1.87759273525725,-1.85097936241481,-1.86934849456797,-1.85711035934366,-1.87373061667475,
			-1.83807642673222,-1.84146984968001,-1.85542091530413,20080.89823704423000,-0.64500465740213,-0.61867298893103,-0.68173241159622,-1.18648741127724,
			-0.77346795519768,-0.38161895743478,-0.38678366425720,-0.38317874445817,-2.29206382118151,-2.19137524805203,-2.27785443629320,-2.25966395554765,
			-2.36717074163473,-3.55537849916466,-2.68475922328713,-2.85966485982210,-2.74061242447461,-2.60056156026921,-2.66745938797056,-2.68626151201305,
			-3.21932200133267,-2.56128970089179,-2.50689898131664,-2.18985331268719,-2.16407672359081,-2.22166487407972,-2.17041405698031,-2.18474753480331,
			-2.15250226743812,-2.17573312862355,-2.16833641968510,-2.19943888191269,-2.14000042939543,-2.15620151401160,-2.13256750960914,19840.71188914010400,
			-0.42698387923362,-0.35975849488983,-0.82309041266894,-1.29676292846107,-0.92290759534829,-0.39397188192860,-0.38636941942224,-0.39147022999142,
			-2.26970908193422,-2.17685430613715,-2.26079916539403,-2.24271039972100,-2.34670281822855,-3.50733140848723,-2.62771582159116,-2.82026932022897,
			-2.72013131519248,-2.57920433305036,-2.63282216171139,-2.66119947976791,-3.19391954590039,-2.53334075609032,-2.46796973365783,-2.14007383500021,
			-2.12471650397049,-2.19354933713997,-2.13654130032160,-2.15577276843614,-2.11273409847871,-2.15159208896917,-2.14056556738835,-2.15521598927344,
			-2.07766517964571,-2.09658227178232,-2.07335957454295,28330.08913184986000,-0.47233627069712,-0.45719424686980,-0.40073087746748,325.38713560290890,
			-0.01766037920414,-0.27234167843045,-0.27032148786460,-0.28111230824086,-2.36533617654922,-2.24332047707050,-2.34780893132602,-2.32060765403106,
			-2.44257226344656,-3.67544413264236,-2.79306539980350,-2.94112821865463,-2.80444390930203,-2.66844956251444,-2.76141094927555,-2.76170219710145,
			-3.30475136579921,-2.63529989136732,-2.59990575821973,-2.26660774649118,-2.26006951293217,-2.29952384509779,-2.24172854050794,-2.25989439930644,
			-2.24298949534544,-2.24548794178436,-2.23235390492587,-2.28470755973807,-2.25019572176689,-2.24735724745277,-2.24485430827113,28323.41271652076500,
			-0.58634209129173,-0.55891108116749,-0.92421768655490,324.74081191358925,-0.67553323317585,-0.36884845354671,-0.38057756580088,-0.37345305562804,
			-2.37076437356517,-2.25563787574136,-2.35506934350316,-2.33460267390675,-2.44600726560863,-3.65799479226701,-2.77023146640792,-2.93066129842234,
			-2.80926939330376,-2.67075329282443,-2.75294010616680,-2.75824077073291,-3.30178436572911,-2.63391839435662,-2.59432145793491,-2.27190048911024,
			-2.25495049562303,-2.29828644290218,-2.24624505880779,-2.26011669006455,-2.24045159137355,-2.24987740008402,-2.24215580276543,-2.28225661578592,
			-2.24008890226276,-2.24999411499590,-2.22978295459135
		});
		return emat;
	}

	private static EnergyMatrix makeExpectedEmatDiscretePROtoPRO(SimpleConfSpace confSpace) {
		EnergyMatrix emat = new EnergyMatrix(confSpace);
		emat.fill(new double[] {
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY
		});
		return emat;
	}

	private static EnergyMatrix makeExpectedEmatDiscreteLYStoVAL(SimpleConfSpace confSpace) {
		EnergyMatrix emat = new EnergyMatrix(confSpace);
		emat.fill(new double[] {
			676192.06292413160000,15700210.44002362200000,947.96476807853800,563058.63400449900000,557.81091092044280,56.31146974464652,125.54663891715616,6.94928834000398,
			7688.92492005508200,14.62258043235105,-18.06128231515508,-17.81097550295133,-12.00693222360090,-17.39842345933905,-17.19537662018083,187.23023627666896,
			-16.04010817202500,-16.25544859314653,-15.48973720840127,-16.49238653156751,-16.08018918212045,237.10972166587527,3336.45825670001570,1339.54650766605030,
			1655.00600504863840,1478.39968741484700,617.69184809913880,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,11.17991322231877,-0.20102714862923,-0.13590859711377,-0.19798909154113,-0.20202557801931,
			-0.17172414296713,0.03224608437460,-0.15643731440255,-0.15850104720200,-0.04811394772883,-0.20028241493973,-0.14177775390480,-0.17676719589748,
			0.00065133230863,-0.22652539993474,-0.21207055040743,-0.20012744691825,-0.17261048893338,-0.18441902015163,-0.19438953043635,-0.18978253112663,
			-0.19025631462282,-0.18955722543661,-0.19512494309969,-0.19469553500860,-0.19563463762425,-0.20111577957099,-0.17355616191361,Double.POSITIVE_INFINITY,
			Double.POSITIVE_INFINITY,-10.82447883853477,-0.23940375591853,-0.15417041270260,-0.22441898801964,-0.24194068032448,-0.19337700602729,-0.00114149801872,
			-0.19702763679167,-0.18957224633576,-0.07498357793164,-0.23955672788405,-0.17948043872501,-0.22936128902470,-0.03179604875354,-0.27269584855873,
			-0.25503453565345,-0.22936857883670,-0.19149202104523,-0.20477715621966,-0.22392861042000,-0.21947155336488,-0.22008576734290,-0.21981969804813,
			-0.22510629963073,-0.22430486344554,-0.22524474677125,-0.23053040484466,-0.19312208324558,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,10.01987919567381,
			-0.25352364943301,-0.15286843580880,-0.23909348833652,-0.25715519204081,-0.20841312350872,-0.01866228851837,-0.22059676038493,-0.21435951497851,
			-0.08515680856555,-0.25292907329032,-0.20526808293505,-0.23979461927714,-0.04622086656177,-0.28479338020391,-0.26700609971830,-0.24077764806474,
			-0.20254088495819,-0.21591682632502,-0.23532410639514,-0.23075080060925,-0.23130630747196,-0.23147834092264,-0.23685556741272,-0.23550500064339,
			-0.23637308736941,-0.24185548094836,-0.20391913397947,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY
		});
		return emat;
	}

	private static EnergyMatrix makeExpectedEmatContinuousALAtoVAL(SimpleConfSpace confSpace) {
		EnergyMatrix emat = new EnergyMatrix(confSpace);
		emat.fill(new double[] {
			-11.18489763171806,-12.90750796191863,-1.27210172099789,-13.62032555012168,-1.40544228323570,-12.41811771090904,-1.44468137674953,-13.88532077176064,
			-1.29877676876898,-12.58059681909518,-1.44496979801170,-12.49358120555075,-1.39453328825675,2.93929493090544,-1.03695453435824,2.79231607127480,
			-1.09476163963637,2.67172360958455,-1.05952279026461,2.94034779784064,-1.07106315842574,2.72603810816966,-1.08343785962412,2.66449260040225,
			-1.04083459419467,-15.42949086105453,-1.16839173114295,-16.11737388387265,-0.88696402213489,-15.65316984762600,-1.09469269353261,-15.92416174272980,
			-1.08393021290948,-15.92593291248491,-1.03165899698177,-15.51703431121919,-1.13802487283863,62148424.45101786400000,-0.25165654256473,-0.35798661660065,
			-0.33600194863105,-0.37098104634105,-0.34927159200795,-0.35799361040597,-0.37031506674894,-0.70947352540963,-0.42578631332648,-0.64553815245269,
			-0.56096729106009,-0.58432890006060,-0.69086960667874,-0.66020812766203,-0.74968102292286,-0.77192179457818,-0.68719897191782,-0.77611385231629,
			-0.72919616253937,7107.13354266541050,-0.27474866796674,-0.42269402961463,-0.39719467808747,-0.43497969456690,-0.41096568787487,-0.42298445421164,
			-0.43382790738811,-0.75757258600901,-0.54884002330061,-0.68420744462024,-0.62134842468614,-0.63638430715266,-0.72926059145973,-0.74325124800099,
			-0.80790812560924,-0.82901011185598,-0.75892231066550,-0.83040372888016,-0.79866341792947,31.01362231529220,-0.23588663831096,-0.42163387962477,
			-0.38918164839753,-0.42064200466929,-0.41102523558739,-0.40240862413654,-0.42625118014841,-0.75867372497870,-0.54637316247299,-0.69302891844917,
			-0.61975412336352,-0.64188915151212,-0.73454518490225,-0.70851982198590,-0.80079985084384,-0.80600832952751,-0.74837280288847,-0.81609132903318,
			-0.76238217360840,28.96549219259054,-0.23741093915591,-0.43337607255304,-0.39776867892259,-0.42623275476746,-0.42493194638681,-0.40953932235471,
			-0.43762954194278,-0.77578795889500,-0.56368961797004,-0.71081533586030,-0.63289036345033,-0.65665953649890,-0.73944998019709,-0.72276304814299,
			-0.81577110223225,-0.82265944807554,-0.76231790838541,-0.83211197752499,-0.77831576573367,3993478.82330598540000,-0.26087022837716,-0.35865260181078,
			-0.32570731225440,-0.36020151688809,-0.34706024957689,-0.34570184004228,-0.36465613434783,-0.71180819892190,-0.47736201036967,-0.64595541065624,
			-0.56763106096602,-0.58514426687115,-0.69108107531993,-0.64541002086150,-0.75030183348236,-0.76137751554332,-0.68741243960208,-0.77094150118246,
			-0.71240547717300,11.63276691450083,-0.21820706207673,-0.46900889420742,-0.43669227196847,-0.47447960332661,-0.45722707692036,-0.45943771542545,
			-0.47686513527855,-0.79477428286187,-0.60061062797369,-0.72157187249577,-0.66657484530169,-0.67906334504179,-0.76472079455987,-0.78353043659531,
			-0.84720502389097,-0.86088391367517,-0.80383023396861,-0.86434385743563,-0.83230534569784,13.25067265395476,-0.21447377162437,-0.44575108077665,
			-0.41167507832382,-0.44979205834076,-0.43395144104242,-0.43104294641125,-0.45354487190876,-0.76836738903812,-0.58080891058564,-0.70559924582417,
			-0.64017721900452,-0.66483359746084,-0.74259802176466,-0.74728688244455,-0.81954600595621,-0.82842078049974,-0.77494081419595,-0.83450584891751,
			-0.79498803478255,-29.83047566613474,-0.18136900488708,-0.46022970590426,-0.42567578572959,-0.48759286269884,-0.43322744071561,-0.47443066083077,
			-0.48135750696116,-0.75355332326140,-0.56227309441488,-0.64426792205485,-0.65924459206314,-0.61326490265307,-0.69512830266941,-0.71148820740947,
			-0.73262699475122,-0.75841350020570,-0.70485753543811,-0.75208010315828,-0.74922952367636,-30.32331879128156,-0.17051332017460,-0.48603434292047,
			-0.44958055280999,-0.52315614329735,-0.44997746679451,-0.47722908151328,-0.51264926752546,-0.75509810660825,-0.52400365495779,-0.59420108993947,
			-0.67011549527079,-0.56833199182899,-0.69538945755559,-0.71658314125058,-0.73367393941417,-0.76194120404793,-0.70735153226587,-0.75384084289501,
			-0.75413395938044,-31.75731135116694,-0.15981158309355,-0.58158852663456,-0.50705681247926,-0.61310269897210,-0.53488034336663,-0.58136129167103,
			-0.61068125996413,-0.82315560883946,-0.66481523048111,-0.72892788385587,-0.74594433702972,-0.71043201883310,-0.76444267948112,-0.79338988352074,
			-0.81045819134719,-0.83144982864154,-0.78816783188036,-0.82537239565211,-0.82523715351457,-30.52103680123505,-0.16186612666336,-0.52093735769447,
			-0.47044900503609,-0.55183273590909,-0.48555328726121,-0.52337501458102,-0.54597123614296,-0.79766476563046,-0.60979175269176,-0.62718272150624,
			-0.71413871892258,-0.65524870419560,-0.73403392035053,-0.76425783960012,-0.77691465792358,-0.80448062978110,-0.75406188085888,-0.79542738869810,
			-0.79867718652496,-30.99361972869080,-0.18917552896825,-0.38234439948142,-0.34470000335767,-0.39008909590730,-0.36877717648191,-0.37450679351683,
			-0.39215719365394,-0.68438710519957,-0.49583987545046,-0.57799791666099,-0.58656501367501,-0.54952519723087,-0.62731409926198,-0.60786074887313,
			-0.64544845220365,-0.65642674022440,-0.61790601744502,-0.65702507900704,-0.64000956844866,-31.82540510553860,-0.19187219202894,-0.37244549655179,
			-0.34256762888181,-0.38067371770748,-0.36091048469391,-0.37021383950981,-0.38075174514425,-0.68233821548372,-0.48782984502396,-0.57485417810292,
			-0.58363894198422,-0.54372077119467,-0.62659189793524,-0.62024300014369,-0.66192599082162,-0.67596001192366,-0.63014828918498,-0.67603677282988,
			-0.65813154830078,-25.40229396662286,-0.20911245349247,-0.32132690465597,-0.28946512371699,-0.33038476542354,-0.30864582359890,-0.31895349506942,
			-0.33182029170821,-0.63114532787149,-0.42888564393214,-0.53344821550855,-0.52310727445432,-0.49582566963296,-0.58208278648449,-0.54933898954001,
			-0.60095283334835,-0.61793899871827,-0.56186006974604,-0.61878361853849,-0.59459971749240,-31.29524762159996,-0.17629829197006,-0.40587804445285,
			-0.37476833924613,-0.41733838993547,-0.39060742426807,-0.40956684028444,-0.41553806737285,-0.70495219058094,-0.50936928889432,-0.58779294057997,
			-0.61542431727788,-0.55678084795213,-0.64293079595028,-0.65859742952145,-0.68479317526829,-0.70183112998534,-0.65935903778333,-0.69872793080335,
			-0.69191739107949,-32.00600692013550,-0.17588591888032,-0.39833110659332,-0.37764596450016,-0.41194384407064,-0.38767640105093,-0.40650205953954,
			-0.40867555444570,-0.70925779931715,-0.49997829442643,-0.58157612632883,-0.61435335537576,-0.54768460861018,-0.64310170397138,-0.65424771259091,
			-0.68098277448994,-0.69935697447584,-0.65372866277405,-0.69577760114476,-0.68925773853942,-30.99983150650667,-0.17297744514910,-0.40429348305383,
			-0.37139643280519,-0.41564620107504,-0.38977773903339,-0.40193765563432,-0.41495458036791,-0.70283929756140,-0.51186288464471,-0.58729523246215,
			-0.61238023207061,-0.55933725365583,-0.63805004382681,-0.65032236411394,-0.66699387968479,-0.69418762921928,-0.64442169812229,-0.69209478463132,
			-0.68326759759871,405.96376369688480,-0.20748719933661,-0.32576904096647,-0.30314933877583,-0.34566038253748,-0.31117049827647,-0.33742364356785,
			-0.34061649840669,-0.63950128923218,-0.42683667284889,-0.53239451188475,-0.53040104022149,-0.49110441952978,-0.58800110795360,-0.57443018523811,
			-0.61331634323136,-0.63855273998901,-0.57555661738753,-0.63511984275893,-0.62087949230993,-31.33926762182741,-0.19743620591625,-0.37876089323576,
			-0.35751578810741,-0.39720017883182,-0.36491515164059,-0.39041270721792,-0.39258463785259,-0.69047206128170,-0.48102405907334,-0.57203379067990,
			-0.59024405379449,-0.53487424442846,-0.63063784921607,-0.63712588413230,-0.66643239706788,-0.68905241875293,-0.63572834298620,-0.68494619258896,
			-0.67654469253747,-34.04736242315079,-0.19531378400383,-0.38545654858649,-0.36055398471969,-0.40396812737864,-0.36986725970784,-0.39542784463104,
			-0.39953230893876,-0.69177841162507,-0.48782615873929,-0.57519229483796,-0.59429220940595,-0.54045981361694,-0.63154744906564,-0.64134951137664,
			-0.66890114087350,-0.69110686688871,-0.63944140349205,-0.68672987360254,-0.67966214743125,-14.51679287725569,-0.24659179625454,-0.47200209389494,
			-0.45275127933756,-0.49712770203636,-0.45443425655466,-0.48967870532044,-0.49094109504323,-0.78614362848486,-0.57135732621308,-0.67738927907258,
			-0.67550144161109,-0.63451501460389,-0.73507542907334,-0.73460929817432,-0.76790219489240,-0.79788084759461,-0.72916403772463,-0.79252009061780,
			-0.78268175581360,-17.03117181333878,-0.25208259598552,-0.48381283730745,-0.46505442788235,-0.50510973201080,-0.46887217389000,-0.49937539754513,
			-0.50003019334443,-0.80139606034280,-0.58242506118933,-0.68665228927492,-0.69022067280958,-0.64506901290559,-0.74897317538118,-0.74664374195229,
			-0.78084841578446,-0.80825403163225,-0.74392214706469,-0.80573064111487,-0.79309078967588,-14.78277173391042,-0.25034562187683,-0.46637918826248,
			-0.43940042035237,-0.48567008082898,-0.45019886525929,-0.47702260788574,-0.48163324612597,-0.77817773044481,-0.56776062694227,-0.67375983082638,
			-0.66645772549623,-0.63322250270570,-0.72843893156030,-0.71816956124460,-0.75742615675765,-0.78207549998687,-0.71953267038247,-0.77884910772660,
			-0.76446373152478,2006.49278406619760,-0.33727571148753,-0.39239827115309,-0.38215884568592,-0.41961720999749,-0.38036856480401,-0.41611425914930,
			-0.41180746016417,-0.73785835963476,-0.49014537821142,-0.63054857850255,-0.60438835865796,-0.57421963848103,-0.69321617680113,-0.65528661306208,
			-0.70753693477624,-0.74665884358856,-0.65227568568762,-0.74130434389560,-0.72102706609635,17.67222370605316,-0.27504838817324,-0.46638830020486,
			-0.45720953233886,-0.49221808005416,-0.45424899531295,-0.48897251348367,-0.48500625539147,-0.78446002360764,-0.56265642967323,-0.66512736982503,
			-0.67746289251648,-0.62886997762797,-0.72962262943964,-0.73713997719744,-0.75879363587432,-0.79188959751145,-0.72927145285751,-0.78578982573760,
			-0.78810007090266,12.93323848766566,-0.26716229179979,-0.45731377646457,-0.44505367716562,-0.48612645396414,-0.44173886123606,-0.48153517772308,
			-0.47886765605125,-0.78154769654909,-0.55301818480333,-0.66831235389778,-0.66438146572854,-0.61962146131047,-0.73042548343677,-0.72552023962113,
			-0.76006454951213,-0.79577570679858,-0.71633570691792,-0.78930128866774,-0.77788782744039,294.63254774185970,-0.29114861469934,-0.41876116877435,
			-0.39602691561598,-0.43447385847697,-0.40913382364127,-0.42779327998627,-0.43247028204216,-0.74594972292520,-0.51985165721079,-0.65281280899764,
			-0.62111949659062,-0.60043890846075,-0.70546953531142,-0.66096892635158,-0.72219842400745,-0.74403414288778,-0.67427769394806,-0.74474863987046,
			-0.71552015229240,466.58259554285300,-0.29358402974773,-0.41787379746772,-0.39482303825208,-0.43381225815279,-0.40578884556258,-0.42495651551131,
			-0.43077431622171,-0.74967552797931,-0.52277410441088,-0.65764666467170,-0.61752095649565,-0.60712428549751,-0.71086067254070,-0.66695188957963,
			-0.72934675523722,-0.75304367469233,-0.68056126379528,-0.75321320458311,-0.72279446889751,210242.45595959094000,-0.30767113779369,-0.40786746386906,
			-0.39202844254251,-0.42874351641299,-0.39609294081194,-0.42351044591286,-0.42386245801039,-0.73964483641413,-0.50095281562767,-0.63453607002328,
			-0.61222382917139,-0.58129221155021,-0.69473687622070,-0.65669668099960,-0.71192863070725,-0.74349568097834,-0.66018828910364,-0.74079445294661,
			-0.71686312271177,-28.17921060598944,-0.22972602306944,-0.59203570796144,-0.57890740987838,-0.63786135877636,-0.56372118794988,-0.63084610630701,
			-0.62507857906585,-0.91031629104669,-0.68937758129805,-0.80647625388854,-0.78941618785864,-0.75967134099074,-0.86398038674139,-0.86501122070675,
			-0.88877371266799,-0.94330561270386,-0.84332002624937,-0.92957261731733,-0.92785323949413,-27.95064833690564,-0.22553803185337,-0.61320666940559,
			-0.59469880075561,-0.66718293772237,-0.57459791748771,-0.65879712696352,-0.65445915718893,-0.91587093840530,-0.69854721393211,-0.81269582550117,
			-0.77362495470731,-0.76843299256167,-0.86835646385179,-0.84105902539050,-0.87063921109535,-0.90626035073727,-0.82593188870808,-0.89324299869893,
			-0.89241672137808,-17.73909844304809,-1.00106323322973,-0.72169198699370,-0.72296068213833,-0.84083629264278,-0.49451264411409,-0.82357765047156,
			-0.81137840759335,-1.03218344664729,-0.81633308702245,-0.96937738939573,-0.88114898399322,-0.91578331757930,-1.01035756435736,-0.82715039953761,
			-0.83102720433507,-1.11940134098660,-0.77779304414167,-0.88963030338891,-1.10072621382546,-30.41892811823278,-0.21762742382318,-0.58329711875390,
			-0.55782228729910,-0.62113895289059,-0.55044153301744,-0.61059143140551,-0.61218286545840,-0.88235171599817,-0.68016208934467,-0.77529313722555,
			-0.77693634561176,-0.73908359109786,-0.82878790065031,-0.84050517404776,-0.86009039052700,-0.95581492013834,-0.82699990214391,-0.94326654329991,
			-0.94492006855992,-30.45521419504089,-0.25136289754176,-0.58816628715119,-0.59831756511878,-0.62263094629563,-0.59436930770682,-0.61313528530562,
			-0.61429979445445,-0.89600956359598,-0.68622161959768,-0.79170108774367,-0.81697774458662,-0.75120390664173,-0.84629988880050,-0.85133733572379,
			-0.87822319188930,-0.91523603086349,-0.87394019959659,-0.90685806432831,-0.90176848600294,67.91390142917194,-0.03336389179825,-0.78300002875006,
			-0.79635716048782,-0.78042735096059,-0.78397395384004,-0.78585550368439,-0.78305031981370,-1.16189944689206,-1.13683733366454,-1.12474306135380,
			-1.15395583443486,-1.12883566117128,-1.13644908613607,-1.00789804609320,-1.00580586763800,-1.00339976734077,-1.00754120134724,-1.00411287783999,
			-1.00577112295775,-1.61078412973303,-1.60763391575158,-1.55663283073910,-1.60009756528723,-1.59916677888007,-1.58450081411276,-1.54105898468508,
			-1.83224892153080,-1.81780536622843,-2.01997205840378,-1.85611437332587,-1.64899407958563,-1.70846375586343,-1.70426311841951,-1.66492843799793,
			-1.68050025457914,-1.63698759193754,-1.71274919902120,-1.71297476056467,-1.68193185997541,-1.47707978372165,-1.48448348973802,-1.46386342640852,
			-1.50128155241910,-1.49593306098728,-1.48760547912361,-1.48361636853447,-1.50216284538221,-1.50184338794290,-1.52983475089816,-1.50287870483952,
			-1.51329673998696,-1.47807706839558,-1.50533279515796,5.33463145681181,-0.03931005664915,-1.04727602272135,-1.00500826616579,-1.04201208110555,
			-1.03740960566862,-1.02301807583928,-1.04765056396195,-1.29453744618942,-1.30105937430404,-1.28332175151429,-1.31584302907084,-1.29042011695892,
			-1.27143807783508,-1.25094576274005,-1.24903569400485,-1.24718809926355,-1.25068784208583,-1.24723753716724,-1.24934593781243,-1.79043566451327,
			-1.79305988108667,-1.73594705662409,-1.78666195921542,-1.77970318707481,-1.76937547606329,-1.71895390590804,-1.98709408546551,-1.97702888418085,
			-2.14829034620775,-2.01040587211183,-1.81302740758845,-1.87353231837455,-1.86994351358943,-1.83451517011330,-1.84381504020696,-1.79963580998639,
			-1.87323253735775,-1.87683282295549,-1.85190398085677,-1.64854621812129,-1.64941013261388,-1.62337676025993,-1.66633155357925,-1.66031415216379,
			-1.65317750201442,-1.64394886176867,-1.66792311217805,-1.66801302569918,-1.68901690495652,-1.67102092128813,-1.67324885384341,-1.64438951091362,
			-1.67490703905393,26.55919340130361,-0.04636165466968,-1.05564771889399,-1.01905955363948,-1.05353319247568,-1.04551461478289,-1.04571513429260,
			-1.05658779970506,-1.24480922296165,-1.24071007150578,-1.22438294466000,-1.25862064582877,-1.22966237801246,-1.22307542488975,-1.28934517580578,
			-1.28742975752912,-1.28597748045408,-1.28900480646252,-1.28594565227356,-1.28792466030451,-1.95527961679578,-1.95258281485282,-1.89396090427835,
			-1.94529871889758,-1.93882733210617,-1.92855776453880,-1.87721173478821,-2.16992423191256,-2.15428417158210,-2.35499070911050,-2.19525280019923,
			-1.98067339181818,-2.04155898098814,-2.03847925959479,-1.99712315180896,-2.01257957462762,-1.96765894908132,-2.04668603904241,-2.04589851764881,
			-2.01489250218969,-1.80549748528581,-1.81262000467956,-1.78690080541679,-1.83053267696284,-1.82441466959772,-1.81638952572741,-1.80788541644545,
			-1.83094898952420,-1.83141347639465,-1.86276432649849,-1.83566792046223,-1.83924510466491,-1.80112102458867,-1.84671875609545
		});
		return emat;
	}

	private static EnergyMatrix makeExpectedEmatContinuousGLYtoGLU(SimpleConfSpace confSpace) {
		EnergyMatrix emat = new EnergyMatrix(confSpace);
		emat.fill(new double[] {
			-10.65318811007568,-10.12769381224124,-0.99421108301495,-10.83084871997190,-1.13436220064404,-9.50314073726129,-1.19995922267993,-11.26849222565858,
			-1.00170935818491,-9.68987428118261,-1.19477543656326,-9.56914985547872,-1.15252505789243,-12.98190832601781,-0.95883511982553,-13.21054897069035,
			-1.02009664235172,-13.22895559713672,-0.98281557931366,-13.05327447872613,-0.99989396240022,-13.21799389217471,-1.00691900547225,-13.20677083877471,
			-0.96294223803743,-12.95150915481017,-1.24360361704375,-13.60144013602408,-0.97169830245614,-13.20686416624782,-1.16604098605920,-13.39584125521438,
			-1.16964620503581,-13.44946072565427,-1.10560873262216,-13.05348859613299,-1.21004903074684,-3.39456151481206,-0.15909213643326,-1.70035870083883,
			-1.79227683574721,-1.72032242563131,-1.74839741630364,-1.77288623690988,-1.68919864239207,-2.07953384223109,-1.80860484848992,-1.43451631895305,
			-2.18943832827726,-1.47283336821920,-1.65142378065644,-1.78546075149148,-1.70661042738182,-1.73094743317737,-1.73615377485467,-1.70207382611278,
			-1.77689386867578,343.98486711647650,-0.17268289498630,-1.49411893564938,-1.59563269820590,-1.52548013820116,-1.54204737448037,-1.58099662308794,
			-1.48913579952996,-1.88732631015057,-1.52092685427047,-1.30782370066440,-1.92974260963630,-1.27967933940378,-1.53240000387884,-1.57364343073510,
			-1.52157533942012,-1.55382971079886,-1.53283101404856,-1.52902447670294,-1.58557779016296,-19.60006734063828,-0.17179055176688,-1.17258434274907,
			-1.24349257309115,-1.18674032442351,-1.21940259615814,-1.22637096718611,-1.16330458405824,-1.45709651159064,-1.12332448880445,-1.15351281319303,
			-1.43498715730476,-1.09825259081964,-1.28994600319667,-1.21555855165286,-1.21233510181948,-1.23390845355675,-1.20477272833156,-1.22087617113532,
			-1.23723218607286,-19.18574732452036,-0.17147493899212,-1.17063225767644,-1.22448524553946,-1.16896367249703,-1.21563014888341,-1.20091786008439,
			-1.15387320903333,-1.43522394019322,-1.12545460263793,-1.15529707177352,-1.42064456975572,-1.10779516088249,-1.27545263859394,-1.19899489235429,
			-1.20435230990065,-1.21602340931129,-1.19679474347075,-1.20863436912872,-1.21714707071964,12.20326797713824,-0.19154769595043,-1.47709642589037,
			-1.52712834876909,-1.45207223481670,-1.52857801618941,-1.48826510922152,-1.44502779572653,-1.73533555600830,-1.43240444505183,-1.42482924949645,
			-1.73740611310166,-1.38417754003475,-1.55173394635510,-1.49978986547905,-1.50396516411312,-1.49385274199992,-1.50869109375961,-1.49073280794827,
			-1.50160582840055,-15.09137857771600,-0.21877273242620,-1.24772311062531,-1.28942540354142,-1.23991101664050,-1.28683075410998,-1.26799652125079,
			-1.22910040819749,-1.50682051417718,-1.18724143624090,-1.24819807011052,-1.48171617145139,-1.19660199032387,-1.36178228535655,-1.26276102012230,
			-1.28506292537347,-1.29173845048248,-1.27151076536890,-1.28875963232930,-1.28748699961489,-13.76379517811174,-0.11917424333239,-1.39560617447465,
			-1.43342304006085,-1.40402564495961,-1.41624169251626,-1.43004931797896,-1.39069441836659,-1.62365438464506,-1.36990898614216,-1.32353465346173,
			-1.64058484100647,-1.29604770252894,-1.44143639939303,-1.42448657126590,-1.39738327278129,-1.40617254986661,-1.41192298698204,-1.39221971571004,
			-1.42518378687154,96.05247496923565,-0.11894356893001,-0.24630284918594,-0.45971424018075,-0.25996468882890,-0.35987178743062,-0.42664102964786,
			-0.17799901250180,-0.50821234862410,-0.27456804893294,-0.36798897704396,-0.37679129113368,-0.27988957409492,-0.48070190021613,-0.37580931504970,
			-0.41809917226501,-0.42569775946106,-0.37823015901590,-0.44008544676743,-0.39902188889028,-1.35572129328627,-1.35339099664136,-1.33797935404155,
			-1.50125591470508,-1.35632911103849,-1.40681901629267,-1.50310616202999,952.62456137483910,-0.22141836090803,10.89619717039654,10.37627780656979,
			11.07794190666963,10.42945075803837,10.76660741789163,11.14345011223206,10.49377435137128,10.88643767764811,10.86697738565151,10.65729236479607,
			10.95091442818097,10.67323084774037,10.71411848929809,10.67999311639162,10.72020511293795,10.70198430529404,10.68918184795531,10.72432321794355,
			-2.26177890017693,-2.24827554669552,-1.99306585385805,-2.17773702613540,-2.13578871902358,-2.27697068646072,-2.68801573108619,679.30697590741480,
			-0.05174847916229,-0.29806248491534,-0.70824412103018,-0.49884344201996,-0.40317290815484,-0.81595404321495,-0.31158198259553,-0.71649522866535,
			-0.24663114597827,-0.44915722937737,-0.40507798064769,-0.24613003788481,-0.80627540181456,-0.51751420717403,-0.48012214764099,-0.57547787809363,
			-0.46703519261623,-0.53711443844106,-0.57728859783185,-0.82139136517046,-0.72234597036381,-0.88744374924347,-1.09991240655333,-0.96609976678267,
			-0.99227724942569,-1.01429242346767,1738.43928583474830,-0.07573189654350,-0.14884405866631,-0.45637341996669,-0.27631924529960,-0.26908006918544,
			-0.51685707781977,-0.13643226351151,-0.45947965904065,-0.20649173494561,-0.33741938901212,-0.30527031151482,-0.22155597558696,-0.48417362115584,
			-0.32542666237300,-0.35189768973042,-0.40081102704260,-0.31617140291846,-0.40636063064631,-0.36686189737593,-1.16341915375973,-1.14512827177334,
			-1.16677453872773,-1.34210825083321,-1.19972211431136,-1.23177564764818,-1.28337457125928,5782.52616105736000,-0.06464112210692,-0.32726696311201,
			-0.50604023390386,-0.36235572212533,-0.42645893413193,-0.46025571561500,-0.30215189773799,-0.50898650928530,-0.32217317550380,-0.36641671089396,
			-0.43015403858025,-0.30963588367364,-0.45238786667948,-0.42197069028581,-0.42365419502435,-0.43029645788634,-0.41213054912196,-0.43445149076007,
			-0.43181197358633,-1.27622811032044,-1.22140891285110,-1.13461522424306,-1.33832911224930,-1.27203849020877,-1.29081835261714,-1.40415333415733,
			-17.04074847997199,-0.31038732166242,716.90484597389140,706.51077587152700,713.52546687930950,714.76248673040200,697.69083757344340,716.79693776521700,
			301.99585915141114,303.57039898324683,303.00886800879510,302.99552089915870,303.58912360174367,301.93088602901410,242.56911192396880,242.45278382992427,
			242.34265926830530,242.56370979404810,242.34215566808854,242.45860861519557,-1.99501137704554,-2.06543458457506,-1.98164272029560,-2.17360009089905,
			-1.99513330473940,-2.18487023667864,-2.39307856214536,-20.50523390579272,-0.10765403113634,-0.33788344152870,-0.97446550706648,-1.52999299319829,
			-0.48745394390690,-2.12982756120620,-0.62802609955081,-1.15823737626720,-0.59926212952920,-1.17173756064742,-0.73308316900754,-0.74678808296524,
			-1.48972845092098,-0.84047767290433,-0.92736331084236,-1.09363805993511,-0.80946820311883,-1.09457051400972,-0.96867283762431,-1.35815823722845,
			-1.37629645207363,-1.43275569107812,-1.60596129897298,-1.44994492674975,-1.50713005950168,-1.53734229093257,11.30201067304476,-0.06259292600592,
			-0.71478258571281,-1.28549903656272,-1.06520075095880,-0.85471262345867,-1.39298577925593,-0.78244388476723,-1.53181903654049,-0.49779920126298,
			-1.10309439285200,-0.75272298573154,-0.59072472098148,-1.86932811437742,-0.99692259458348,-0.91497642554833,-1.05057782806965,-0.90994432195492,
			-0.97646057888713,-1.07784026187911,-0.92560951136946,-0.82357990598751,-1.08980183718311,-1.30093862598072,-1.15656096448173,-1.20290733336528,
			-1.19435811104243
		});
		return emat;
	}

	private static EnergyMatrix makeExpectedEmatContinuousASPtoPHE(SimpleConfSpace confSpace) {
		EnergyMatrix emat = new EnergyMatrix(confSpace);
		emat.fill(new double[] {
			73.67471122657665,304.41789129678233,-20.33174505008208,-21.51027989888533,-16.69487733379109,6.84281318533537,-1.24857666150071,-1.27595763421853,
			-1.17686196694510,-1.24958347286303,-1.18928876969133,-12.48563569307660,-1.07651408950785,-1.10307230117628,-1.10601941408836,-1.13837203711853,
			-1.08587106974603,-12.79329722361563,-1.11354978451486,-1.15619264297584,-1.23073968675085,-1.27335096216550,-1.18774219520696,7484.47146699407500,
			-0.25643911310533,-0.25525122419378,-0.22076835386308,-0.25193328464788,-0.25888806356821,-1.33605530528500,-1.45628398864490,-1.32866622544925,
			20.50100502536608,-0.17509475687465,-0.16916829681402,-0.15888657479308,-0.17528907132589,-0.19624623800152,-1.33090432933384,-1.44253579159920,
			-1.31394731086840,-11.58169663017034,-0.22309469327526,-0.21165876953852,-0.14929048544461,-0.19335460556112,-0.22991896648933,-1.22674508746964,
			-1.29113306586144,-1.22069190639851,-11.45894095309919,-0.16421956690412,-0.16351660119268,-0.07360722266950,-0.09378839062327,-0.17360458435069,
			-1.21227067737454,-1.27430338814733,-1.20436004422038,-11.75767291725103,-0.18697822833385,-0.18668443486435,-0.15598265432195,-0.16877476010685,
			-0.19330111122393,-1.19238398679101,-1.26186835117242,-1.17481224635777,-11.48310558991511,-0.16945385153155,-0.16604735999378,-0.15192685282265,
			-0.16344208345288,-0.18644399555040,-1.19629634359098,-1.25090916623615,-1.17882469469725,-11.83946782931571,-0.22024618900752,-0.22287249819232,
			-0.19116682820417,-0.22011288821670,-0.22443306865263,-1.19738900915297,-1.25321136944888,-1.18102047863286,27.20162948808694,-0.19910518646980,
			-0.18916334425615,-0.19045318320507,-0.20704890944912,-0.23631041513598,-1.32862371090187,-1.38740837560811,-1.31142588669794,31.01971832373620,
			-0.27715968151925,-0.27953703467334,-0.27463466026326,-0.31576222140836,-0.29355589699322,-1.32759654660152,-1.38740065465701,-1.31625623313156,
			2759.36860992814500,-0.25678041798797,-0.26941183283525,-0.23377944091367,-0.25905635978890,-0.27763249839254,-1.32734820341794,-1.40744547868900,
			-1.32103666493814,-12.54386278830734,-0.73506121287442,-0.79460727658204,-1.27413515152292,-1.37301134753400,-0.65775934344140,-1.73311101378284,
			-1.79214412829377,-1.71021167103008,-12.03968212468207,-0.40693843055700,-0.40656807749194,-0.44254826291097,-0.64309151236566,-0.36022926802358,
			-1.46288511595852,-1.53254038780492,-1.46287130862178,-11.71804193696711,-0.18790427290744,-0.19444287834147,-0.05211621173207,-0.06975788822153,
			-0.21447225856946,-1.45410326979650,-1.52327770552140,-1.44601342454443,953369.55727655720000,-0.10491716454632,-0.09964351319178,-0.12046705085529,
			-0.10386444623351,-0.10143754464286,-0.91979564240654,-0.64468407084913,-0.52504161848315,-1.64841842786044,-1.67792670078924,-2.02152086028439,
			-2.18394665185971,-1.77808954834297,-1.77592365424465,-1.76979195156940,-1.57166208588317,-1.55334261983044,-1.55601149204512,-1.80137012925978,
			-1.65670793374392,-1.74639149569868,2669430.04963119330000,-0.03754071236106,-0.03695169556130,-0.03401359420001,-0.04291631324134,-0.03343714202012,
			-0.85984107884550,-0.77385403228769,-0.62439445578135,-1.25957466524995,-1.28901742048232,-1.57350792271769,-1.73293363129997,-1.37469905301599,
			-1.36858127166733,-1.35732748186094,-1.16177793321846,-1.14754663566324,-1.15222281719547,-1.25371134770579,-1.22011811700183,-1.28778634226480,
			49.54898271514714,-0.16713037704362,-0.16896455615686,-0.17252888539954,-0.16454754858090,-0.17316865990533,818622.25427826110000,360.32312158103840,
			178.18919317542650,-1.92378362292525,-1.95243327565212,-2.22693576931846,-2.37803741462926,-2.01624714242141,-2.00144973716891,-1.99321061134741,
			-1.79878129569502,-1.78239769632964,-1.79364336890681,-1.95586914457321,-1.87406601776336,-1.94569878300913,13.21666359124336,-0.12681618926077,
			-0.12369146757690,-0.07733944368357,-0.08846351662022,-0.10154850896519,10.22496544076035,-1.16049723009893,-2.36897372402950,-1.92081568485162,
			-1.94835693011472,-2.19439069216680,-2.34446612096002,-1.99772972717080,-1.98298348213206,-1.97559195965150,-1.77383877496044,-1.75893992361447,
			-1.76733538229770,-1.89923383578647,-1.83305108003509,-1.90255836149893
		});
		return emat;
	}

	private static EnergyMatrix makeExpectedEmatContinuousGLNtoTRP(SimpleConfSpace confSpace) {
		EnergyMatrix emat = new EnergyMatrix(confSpace);
		emat.fill(new double[] {
			177.57340659124046,134471.85107671420000,4.63442862911581,6.29072086080121,-22.61413507928436,-10.84843156055110,-23.02095476159237,577.40115778992610,
			2540.14364337331240,226.53747523811182,-1.73148062961188,-1.78492266651365,-2.06197795905233,-2.29664165953465,-1.84400414022193,-1.61460023480157,
			-1.61788114012730,-1.66567611919091,-1.68919929150482,1735.03248191992220,-1.72170572814736,-1.77517805205693,-2.04012481306383,-2.27455739007734,
			-1.82875143435611,-1.60566416544616,-1.60822828998294,-1.65579615587933,-1.67829414820366,27.97255785520537,-1.43642208832056,-1.47715706521578,
			-1.75136017973062,-1.99696341429978,-1.55323913771031,-1.32250888622831,-1.32877230832403,-1.36977403310350,-1.39445509099936,18.97876924287457,
			-1.38126573617972,-1.42009669035120,-1.69935875183756,-1.94118860534126,-1.49943941912707,-1.27429426127954,-1.28052301888197,-1.32175540161219,
			-1.34588423696221,45.91727856361817,-1.95208685644148,-2.01703295984777,-2.24376929122270,-2.48507348244152,-2.04768096775927,-1.83411466705572,
			-1.83116554298576,-1.87975501838992,-1.90394989814735,18.43696674428394,-1.60088621476196,-1.65781430375266,-1.88214718172778,-2.12632428657886,
			-1.69812805265982,-1.48058639943180,-1.48224171296677,-1.52198688432552,-1.54446522952452,85.10256095104619,-1.79831190060252,-1.87296047451883,
			-2.05873009890847,-2.30693129663475,-1.87697158800205,-1.65159132406034,-1.65305909309629,-1.69173859700122,-1.71503791439334,72773.41275980925000,
			-0.23796197752616,-0.23287659233401,-0.28601663589033,-0.37196709543522,-0.26712341775637,-0.21323920172039,-0.22089785076419,-0.22812954751645,
			-0.21550584635707,-1.99896921148996,-1.66073352629438,-2.47411994543748,-2.09936971325962,-1.92929375045005,-1.93591339750748,-1.98273825043779,
			-9.90413410987563,-0.22665677982623,-0.21748605593703,-0.29065091353618,-0.37661703646208,-0.25845667009781,-0.20091408318370,-0.21088770109140,
			-0.22015567263357,-0.20382872572398,-1.74543478854085,-1.44846212767835,-1.84793941569387,-1.71743790379523,-1.56887192492157,-1.57025020110408,
			-1.58975832219483,-8.22883978367654,-0.21947892541240,-0.21069580156749,-0.27495635413650,-0.35822130568980,-0.25031059146737,-0.19504566439685,
			-0.20458197957799,-0.21384199569631,-0.19754291438794,-1.71409073245872,-1.41157364135574,-1.82463082312750,-1.69148208637083,-1.54186310718574,
			-1.54385144847911,-1.56277714153380,-4.56353924380495,-0.45632664219056,-0.44931898540030,-0.62112217315874,-0.67888271307277,-0.58013546315845,
			-0.36693414335196,-0.38189040695883,-0.39185653770339,-0.37598638860466,-2.09345548972827,-1.79429477986842,-2.27809704512686,-2.07993921537866,
			-1.93872957500156,-1.93591828944880,-1.98772356147998,64.18046394155401,-0.34092175426898,-0.32957383673013,-0.50322883226482,-0.56598711736936,
			-0.45891521668900,-0.28508484497802,-0.29985509479253,-0.30925411268219,-0.29380547557482,-1.98346980870280,-1.68667743054215,-2.13566054833797,
			-1.96635861749610,-1.82109703399738,-1.82029922576727,-1.85651856337231,12939422.18061680500000,-0.05241800307547,-0.05424118031881,-0.27145290669803,
			-0.08991796976603,-0.06326027933039,-0.04361924894597,-0.04495332263956,-0.13638725226116,-0.10818366196111,-1.39892708421780,-1.71612896438624,
			-1.56695923593467,-1.60036330636920,-1.31837984179155,-1.32119460799672,-1.26037522064817,-1.80494337425957,-1.86824443343259,-1.91869405795591,
			-1.53903594075724,-1.56904739256933,8154036.28127955800000,-0.49564097153822,-0.49839871270690,-0.93076057585874,-1.30330265845794,-0.71839968951646,
			-0.51512146199785,-0.52070269545855,-0.88367230489985,-0.75817455949238,115.82438390192243,117.08181778976359,117.32271715572105,117.27162279280262,
			113.23164616697545,113.24409479584345,113.33219906516340,-2.47210269838287,-2.56824216474375,-2.58922150636477,-2.31483833563023,-2.34037239130290,
			3676685.75050887140000,-0.04439661903009,-0.04426618323909,-0.07893935158669,-0.06207622050215,-0.05011492224151,-0.04208822974122,-0.04349554113920,
			-0.05165534144104,-0.04810956899773,-1.27293932055820,-1.61143973266361,-1.32311932680177,-1.59314577004355,-1.14072769505537,-1.16756343575941,
			-1.06339269431472,-1.68377594324966,-1.65522793214929,-1.74888809537476,-1.34355872743029,-1.36434724776280,209466.89450263456000,-0.04032760558053,
			-0.04106389562637,-0.07818682480995,-0.05945728262610,-0.04897393601480,-0.04063450228031,-0.04184939016857,-0.04929677388745,-0.04668666207533,
			-1.15878633507302,-1.23797042989848,-1.04456244866005,-1.12129227846293,-0.89826652625369,-0.90780227840816,-0.84221428565231,-1.49428069635847,
			-1.50558948582904,-1.58753824379935,-1.19854969603765,-1.21947275307952,7.25470972227755,-0.85704839775704,-0.87983501109251,-1.28365927682785,
			-1.21551362820866,-0.94953097376180,-0.87005345104404,-0.85598963818786,-1.00070357778238,-0.96570501362316,35937472.44431083000000,8123489.89920096900000,
			16816044.99891694600000,16816141.50884534000000,40994.99531124058000,41132.73057589739000,17266.66735344366000,-2.69061374656468,-2.64172144100552,-2.70528356264930,
			-2.41843977560968,-2.41025922491052,35406.68424495563500,-0.16390354962313,-0.17619575574242,-0.20116908504564,-0.23117862507217,-0.19178450023278,
			-0.17826786001688,-0.17463787077162,-0.18163208934786,-0.18299752043232,5728.83980953559800,237.63429587358607,6.14298081289999,12.31298626455284,
			-0.27157311218203,-0.46057689522699,-1.86765839335083,-2.00277309969946,-1.90877723778094,-2.00038395353469,-1.60441849860829,-1.62369968371626,
			129590.14644873452000,-0.06937137571418,-0.08635182466443,-0.20671920221969,-0.09844227432243,-0.07598235830713,-0.07385853431967,-0.06224292032155,
			-0.10765267065099,-0.10341469665691,8605072.28797081500000,3975.54473605917560,275.00459891535223,292.07310035831455,0.94108942874926,0.76141753929577,
			0.56096466241219,-2.23161636474767,-2.11788310427033,-2.20367982301416,-1.81315754851039,-1.83187906926435
		});
		return emat;
	}

	private static EnergyMatrix makeExpectedEmatContinuousHIStoTYR(SimpleConfSpace confSpace) {
		EnergyMatrix emat = new EnergyMatrix(confSpace);
		emat.fill(new double[] {
			3197340.10354795600000,2117679.56656603600000,-16.23008537837580,-15.89083323199995,-16.09451079938304,-12.62833842186817,-11.60039542698798,-17.13266735063633,
			40002.96062105375000,-0.15027904734376,-0.14896658929442,-0.41575713715481,-0.30308000860096,-0.28364364975485,-0.13036244365056,-0.13794554213685,
			-0.14444379373526,849302.76923118360000,-0.11747683422438,-0.09942752322716,-0.25624769749975,-0.18559506574622,-0.11789904259453,-0.03140564027830,
			-0.05789697594651,-0.06957133970572,21243.13300915922800,-0.20130839690030,-0.19786835116290,-0.40553578845066,-0.31365175596485,-0.25943097655294,
			-0.12136387143044,-0.14272838752891,-0.15300822355007,14071.64631853276300,0.02870830906764,0.03829927458774,-0.11399502485577,-0.04395737052129,
			0.04351822455214,0.14619944158515,0.12083644952690,0.09878859590248,10517.51420879847700,-0.22095368725646,-0.21242700342105,-0.45046235815431,
			-0.34290144068038,-0.31995299906092,-0.18401544815644,-0.19678423149622,-0.20641747166955,487.68856746415200,-0.13448532257889,-0.05100686651711,
			-0.39195047710177,-0.23386675607557,-0.29481832122245,-0.23671085249400,-0.24071551594477,-0.24780711821428,127.00147820495107,-0.10039116855579,
			-0.02978039512946,-0.46408653488714,-0.31398885958457,-0.36525561474774,-0.27533429333540,-0.27380512218734,-0.27788874075933,64.52707950872144,
			-0.16641784684596,-0.10216395110119,-0.46434299728273,-0.31924209257967,-0.36545107709601,-0.26502049362698,-0.26668281588875,-0.27100286290779,
			1791.68139274036750,-0.07202190836883,-0.04349555870274,-0.27154505772563,-0.15830121938775,-0.16426486072262,-0.06774342705567,-0.07955396677583,
			-0.08873066335760,961.59029165379170,-0.11272333681932,-0.09819921286589,-0.34339910351940,-0.23499129468319,-0.23347755168135,-0.10170700641849,
			-0.10836809864209,-0.11637966670055,585.26191551942080,-0.09868075068339,-0.06749562173763,-0.36114000145283,-0.23675251324006,-0.26191621575023,
			-0.14903918955479,-0.15139999257560,-0.15621433417524,10741921.50859303400000,-0.16265630644851,-0.14644052871150,-0.42129917399976,-0.30688895149172,
			-0.32018820688836,-0.18385622434866,-0.18626743142838,-0.19175979088721,653.32056492636210,-0.11973136393232,-0.09000857351548,-0.36164915120775,
			-0.24113137120977,-0.26090020698676,-0.15574041532762,-0.16171039868567,-0.16881305396558,152.68700691747880,-0.08069431108154,-0.06935448450744,
			-0.35821130667438,-0.24226143314650,-0.25476767439106,-0.11033628504732,-0.11204490179212,-0.11497061772691,64.61254789274636,-0.01746467108617,
			0.01305908197645,-0.31189433486690,-0.18736094100854,-0.20923664751235,-0.09266474037166,-0.09498734928341,-0.09775218939019,1248.60584225509110,
			0.26309071228017,0.18040987053287,-0.34269839543373,-0.21720602900944,-0.23605811806531,-0.02552640847052,-0.02156256237791,0.04295066441106,
			84.50788234805697,-0.28110689155941,-0.25960203471027,-0.62704845206596,-0.49164074573637,-0.51857027133163,-0.38980187208322,-0.39253987483695,
			-0.38779149661493,399.04224928596070,-0.41887220896841,-0.40717852009455,-0.71655593229824,-0.59423077830759,-0.60872136634437,-0.45809187145260,
			-0.46012788683229,-0.45872690038653,869.45580520530550,-0.07039938400115,-0.10911371792296,-0.47242111016970,-0.35855061260985,-0.35795656856758,
			-0.17258407402144,-0.16168406077697,-0.14553922471570,8680.96631882634000,-0.37226265763198,-0.38762233833374,-0.70781124301737,-0.59026426640713,
			-0.59716982804256,-0.41990047755470,-0.41916457876545,-0.41004262906252,416.76365560615940,-0.21084920044807,-0.20776090909065,-0.60931779774569,
			-0.48241941989032,-0.50296022465614,-0.35332496352828,-0.35117389964238,-0.34019128221940,5553.30076679157600,-0.37768397192679,-0.37701344289669,
			-0.63293277851408,-0.52421678353667,-0.50891763011571,-0.35315934417220,-0.35847796157577,-0.36346658305450,4239.88996511166500,-0.07708192157942,
			-0.09783048560435,-0.36540336728960,-0.26684572460997,-0.23416971388639,-0.04410024198259,-0.03821865236607,-0.04449213386289,227.52456916709957,
			-0.44560209437082,-0.43958227542664,-0.87749720218644,-0.77007076595549,-0.81444271779263,-0.65862628443254,-0.65572587706931,-0.63753149247227,
			2117.41545016493700,-0.37484937710645,-0.32955995813669,-0.98915664903652,-0.81824617004860,-0.87618466482426,-0.76636815048643,-0.76060891389483,
			-0.73518386716247,6459.87061034363500,-0.16188676767423,-0.18237554233962,-0.96027099309801,-0.77413560530111,-0.84903580461814,-0.71406789342640,
			-0.70889777446194,-0.64246984660816,585.04703380901420,-0.49672312963231,-0.43088366951466,-0.93500058809676,-0.78142335869900,-0.83979245202528,
			-0.73403852991092,-0.73326029825748,-0.72529908009203,9708840.02812727400000,-0.28497684881626,-0.26028392399022,-0.34924251562113,-0.31411101018876,
			-0.29999052954190,-0.25633214375694,-0.25619283157428,-0.25677768031279,-2.57151395209234,-2.38682213854355,-2.55277900042457,-2.51993238747665,
			-2.67221206397275,-5.94847714044562,-4.89079270121560,-5.31799357580086,-3.65874732807777,-3.06836432630308,-3.49977632363723,-3.22968391157540,
			-3.98229191797389,-3.06680989036455,-3.06082470240751,-2.45411303799038,-2.44969549794470,-2.48867184037027,-2.41578795831420,-2.42672826952217,
			-2.34569171827339,-2.43056889772557,-2.39594317266492,-2.53251964586002,-2.46386148239683,-2.50786333889139,-2.85842151207450,9708839.97596801000000,
			-0.24930030330041,-0.21442146936367,-0.33613519067507,-0.28940382520183,-0.28684482998134,-0.24953498553624,-0.24848463937724,-0.24908882952667,
			-2.53657033465649,-2.36567433325370,-2.53256241285233,-2.49185026367103,-2.65343238196105,-5.94455739153714,-4.55506283798808,-5.30814992525088,
			-3.66117051697295,-3.04293949755129,-3.48965734148643,-3.21495947634203,-4.00506690122834,-2.99664573406308,-2.92531609306078,-2.37521860855299,
			-2.32816207196925,-2.43191326887783,-2.34551664287313,-2.34913884240461,-2.25527742119047,-2.38583663495011,-2.34911140500112,-2.41611305026845,
			-2.26523269749831,-2.38778818988394,-2.51657665615049,15.33152091778588,-0.19789624656049,-0.18227262050743,-0.24886702830496,-0.23640428926752,
			-0.20663624706056,-0.18609299252264,-0.18582454403767,-0.18598304524378,-2.07920345632733,-1.97887082251176,-2.07525338626364,-2.04182404377638,
			-2.22043179873147,-3.51040708327393,-2.84306752651419,-3.06536704101992,-2.92866310420566,-2.36661530504957,-2.66246583434804,-2.48928496070086,
			-3.21613411187947,-2.33228480905395,-2.26672048348632,-1.90016807783483,-1.92791869990659,-1.94289582846529,-1.88795724848693,-1.89300799697687,
			-1.87629921234783,-1.90731274448570,-1.88776083585426,-1.91336020881273,-1.88518568624712,-1.88277792257523,-1.96185217126705,15.51193325835155,
			-0.21266747161994,-0.19595258213334,-0.28831032188127,-0.27320141212948,-0.24807833780523,-0.20163957038029,-0.20064376408089,-0.20103350154286,
			-2.07128759042213,-1.97166060670631,-2.06528723242161,-2.03830742068016,-2.20669061420135,-3.48623039985628,-2.83222832943706,-3.03477119516787,
			-2.90877554324120,-2.35365991520484,-2.64572597022633,-2.47169191128015,-3.19592563316795,-2.31889863967225,-2.25285134296486,-1.89779136259130,
			-1.92061469790062,-1.93302255743869,-1.88338849749766,-1.88758172458786,-1.87032663659818,-1.90037021542764,-1.88380368720898,-1.90688058607476,
			-1.87862448784138,-1.87924699853859,-1.95233465500523,12955.21457890626900,-0.97014433090786,-1.05726441481465,-1.14288988327196,-1.54565758877833,
			-1.24321593062513,-0.42777027201983,-0.43453042487975,-0.43042526359838,-2.44172755944792,-2.33627908313670,-2.42748124136601,-2.40584043712042,
			-2.57379586844560,-3.91996405661129,-3.23146996575810,-3.52397907730969,-3.32256389669307,-2.74006917700258,-3.02629353325379,-2.86242631449545,
			-3.61131734837668,-2.70729180159588,-2.66972320853970,-2.29741062241889,-2.28168485504555,-2.33802877139813,-2.27528915540852,-2.27299983047863,
			-2.24462882835187,-2.28656147185156,-2.27497907238195,-2.31147739641676,-2.25868392809732,-2.29086844021976,-2.31061250388958,12726.98665331049800,
			-0.66830683926644,-0.68002190356373,-1.16427281205876,-1.59484945335198,-1.28242168354649,-0.45330178236215,-0.45078281020363,-0.45156456348197,
			-2.41334992207369,-2.31486597344924,-2.40573658143527,-2.38306023622900,-2.54732578096723,-3.86393723343796,-3.16244581466353,-3.47714337218968,
			-3.29345593016348,-2.71249795830339,-2.98582937843348,-2.83250317391383,-3.58162646936067,-2.67176181549835,-2.62043834527817,-2.24961422844599,
			-2.23810972856757,-2.30297859939536,-2.23312336252905,-2.23683578153577,-2.19971816669505,-2.25555306115134,-2.24043333587579,-2.25987456613210,
			-2.18594994376703,-2.21942044512389,-2.24236313196351,20155.74459086458200,-0.71930282526575,-0.74180418917093,-0.60087255120431,-0.49506784428472,
			-0.57508427876885,-0.31172612102887,-0.30162827437576,-0.31742084862704,-2.52952474302727,-2.41188519295528,-2.51071617240173,-2.47887392267361,
			-2.66507438267736,-4.06546325226546,-3.37384562521699,-3.59598929669181,-3.40251318318105,-2.81821418329287,-3.13694585006276,-2.94863007123068,
			-3.70174286494997,-2.80414579013489,-2.79343136240778,-2.39060937604274,-2.39575773712203,-2.43588779925800,-2.35428141459041,-2.35930816816384,
			-2.35611429086475,-2.36754409308584,-2.34837682387359,-2.42472437534689,-2.38138209051801,-2.37941248218858,-2.44035143016343,20150.86572712469600,
			-0.94626034440348,-0.96472868926298,-1.01068617887960,-0.86520163359488,-1.20950672121711,-0.41515621337141,-0.43041427686032,-0.42321008132429,
			-2.54093528451763,-2.42661535128020,-2.52264643213151,-2.49741502456063,-2.67387281306223,-4.05656806530119,-3.35885323953920,-3.59393471071648,
			-3.41102265701975,-2.82602424343043,-3.13626392966317,-2.95210336396509,-3.70709330923684,-2.80858229966853,-2.79441244573956,-2.40110928008188,
			-2.39805138406852,-2.41627876560754,-2.36817060681228,-2.36817371807387,-2.36294027692132,-2.37816103005316,-2.36515820985812,-2.43275193787286,
			-2.38238366937812,-2.39385512418768,-2.43788407856936
		});
		return emat;
	}

	private static EnergyMatrix makeExpectedEmatContinuousPROtoPRO(SimpleConfSpace confSpace) {
		EnergyMatrix emat = new EnergyMatrix(confSpace);
		emat.fill(new double[] {
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY
		});
		return emat;
	}

	private static EnergyMatrix makeExpectedEmatContinuousLYStoVAL(SimpleConfSpace confSpace) {
		EnergyMatrix emat = new EnergyMatrix(confSpace);
		emat.fill(new double[] {
			124.53880777042288,8831.33079682992900,641.30851543256600,18251.81287639589600,187.46241325799900,-8.57814966644351,-12.13052703844445,-16.85988314345488,
			44.61944589778682,-15.76762933481837,-18.42419700800705,-18.10717176304254,-17.37632586018095,-17.99897901443023,-18.08285214089562,-16.25464446937923,
			-17.54100490110568,-17.54728486678962,-17.12966144418744,-17.76402038098534,-17.59087523610395,-7.67571439008538,5.68627518156918,37.63144759670099,
			45.12665147556242,179.04299480124223,-9.24206373735060,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,
			Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,9.87205874706383,-0.20820080570794,-0.19237469723874,-0.20529496297847,-0.20891233997450,
			-0.20071466180557,-0.04352280154023,-0.20644354137933,-0.19939741550235,-0.12136461061388,-0.23341236956944,-0.19115380754664,-0.19828542562260,
			-0.05903492394588,-0.23724789655341,-0.22334674407970,-0.20473981159754,-0.19096001862852,-0.18819412047003,-0.19790334924528,-0.19370038702355,
			-0.19483080150848,-0.19488555173828,-0.19989371209537,-0.19985668004862,-0.20143983850768,-0.20681700817341,-0.19429241316838,Double.POSITIVE_INFINITY,
			Double.POSITIVE_INFINITY,-10.84060695172120,-0.25448136851926,-0.20900608487798,-0.25075166480511,-0.25320604093374,-0.23268357814579,-0.07055998244866,
			-0.24919404435436,-0.24062128260433,-0.15707803222568,-0.28402458646562,-0.23695961354033,-0.24985474912812,-0.08210049776154,-0.29334073979616,
			-0.27611753928147,-0.23713865908472,-0.22306411276013,-0.22315778252517,-0.23010430781482,-0.22602971238047,-0.22698745270513,-0.22752849365014,
			-0.23262796499781,-0.23248642522642,-0.23409453748927,-0.23929769268059,-0.22755791991492,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,9.98391851097557,
			-0.26392916452715,-0.21837988110676,-0.24925006216693,-0.26266226266273,-0.24174594967335,-0.08457406130516,-0.25311589496781,-0.26100460996757,
			-0.16332016313593,-0.29156812575104,-0.24160221862395,-0.25589156401061,-0.09177658226242,-0.29977128571652,-0.28252981732235,-0.24439441557114,
			-0.22977723606829,-0.23016756630020,-0.23745127110831,-0.23320597950910,-0.23391662537847,-0.23492768235748,-0.24020199940448,-0.23949675932628,
			-0.24106159740956,-0.24647516418029,-0.23427521427922,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY
		});
		return emat;
	}
}
