#' customized discrete colors
#'
#' @export
discrete_palette <- list(
 ##ref: http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
     # seurat = hcl( h = seq(15, 375, length = n+1), l = 65, c = 100),
 singler =  c(
     "#7FC97F","#BEAED4","#FDC086","#386CB0","#F0027F","#BF5B17",
     "#666666","#1B9E77","#7570B3","#66A61E", "#E6AB02","#A6761D",
     "#A6CEE3","#B2DF8A","#FB9A99","#E31A1C","#FF7F00","#6A3D9A",
     "#8DA0CB","#4DAF4A","#984EA3","#c6c386","#999999","#66C2A5",
     "#FC8D62","#A6D854","#FFD92F","#BEBADA","#FB8072","#80B1D3",
     "#FDB462","#BC80BD","#B3B3B3","#33A02C","#B3DE69","#4038b0",
     "#ee7576","#e94749","#E78AC3","#ff0000","#A65628","#d80172",
     "#F781BF","#D95F02","#E7298A","#1F78B4","#FDBF6F","#CAB2D6",
     "#B15928","#FBB4AE", "#B3CDE3",'#0173b2','#de8f05','#029e73',
     '#d55e00','#cc78bc','#ca9161','#fbafe4','#949494','#ece133',
     '#56b4e9',"#00AFBB","#E7B800","#FC4E07","#FFDB6D","#C4961A",
     "#F4EDCA","#D16103","#C3D7A4","#52854C","#4E84C4","#293352"),
 blindless = c(
     "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F",
     "#BF5B17", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02",
     "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
     "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928",
     "#FBB4AE", "#B3CDE3", "#BC80BD", "#CCEBC5", "#DECBE4", "#FED9A6",
     "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC",
     "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC",
     "#E41A1C", "#377EB8", "#984EA3", "#FFFF33", "#A65628", "#F781BF",
     "#999999", "#FFED6F", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
     "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3",
     "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
     "#D9D9D9","#666666"),
 customecol2 = c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b",
    "#666666","#1b9e77","#d95f02","#7570b3","#d01b2a","#43acde",
    "#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff",
    "#ecb9e5","#813139","#743fd2","#434b7e","#e6908e","#214a00",
    "#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7",
    "#006874","#d2ad2c","#b5d7a5","#9e8442","#4e1737","#e482a7",
    "#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99",
    "#3fdacb","#bf5b17"),
 monocle = c(
   "#016edb", "#09658A", "#B4246E", "#FF9A21", "#FBC83C", "#ff42ef",
   "#216407", "#ff0000", "#483F84", "#995432", "#9834E7", "#35B5E2",
   "#A71206", "#717E8D", "#7762F0", "#C39CFB", "#FFFA2C", "#00FF00",
   "#CAC379", "#8B8A1D", "#FF8380", "#00FFFF", "#C0C0C0", "#0118FA",
   "#000000"
 ),
 col50 =  c(
      "#982f29", "#5ddb53", "#8b35d6", "#a9e047", "#4836be",
      "#e0dc33", "#d248d5", "#61a338", "#9765e5", "#69df96",
      "#7f3095", "#d0d56a", "#371c6b", "#cfa738", "#5066d1",
      "#e08930", "#6a8bd3", "#da4f1e", "#83e6d6", "#df4341",
      "#6ebad4", "#e34c75", "#50975f", "#d548a4", "#badb97",
      "#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344",
      "#d0cdb8", "#421b28", "#5eae99", "#a03259", "#406024",
      "#e598d7", "#343b20", "#bbb5d9", "#975223", "#576e8b",
      "#d97f5e", "#253e44", "#de959b", "#417265", "#712b5b",
      "#8c6d30", "#a56c95", "#5f3121", "#8f846e", "#8f5b5c"),
 ditto = c(
     "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
     "#D55E00", "#CC79A7", "#666666", "#AD7700", "#1C91D4",
     "#007756", "#D5C711", "#005685", "#A04700", "#B14380",
     "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
     "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C"),
 cold = c(
     "#1660A7","#FF6A00","#219418","#CD0C18","#814BB2","#794339",
     "#5ddb53","#5066d1","#FF0000","#F0E442","#4836be","#EE6677",
     "#DC59B6","#CC79A7","#982f29","#AFB400","#11B3C6","#292344",
     "#E69F00","#AA3377","#009E73","#576e8b","#0072B2","#D55E00",
     "#00AFBB","#00FFFF","#4477AA","#228833","#CCBB44","#66CCEE",
      "#56B4E9","#BBBBBB"),
 Buen = c(
     "#00441B", "#46A040", "#00AF99", "#FFC179", "#98D9E9", "#F6313E",
     "#FFA300", "#C390D4", "#FF5A00", "#DBA13A", "#7D7D7D", "#4B4B4B",
     "#8F1336", "#0081C9", "#001588", "#490C65", "#F26C64"),
 pal27 = c("#E31A1C","#3182BD","#FEB24C","#CE1256","#31A354","#91003F",
           "#756BB1","#C51B8A","#DE2D26","#7FCDBB","#1C9099","#FC9272",
           "#ADDD8E","#FA9FB5","#43A2CA","#2CA25F","#8856A7","#E34A33",
           "#023858","#636363","#FEC44F","#253494","#A6BDDB","#C994C7",
           "#D95F0E","#2C7FB8","#FDBB84"),
 pal19 = c("#FF6F00FF","#225EA8","#FED439FF","#197EC0FF","#C71000FF",
           "#FD8CC1FF","#370335FF","#008EA0FF","#8A4198FF","#FC4E2A",
           "#709AE1FF","#71D0F5FF","#91331FFF","#C80813FF","#46732EFF",
           "#FF6348FF","#D5E4A2FF","#ADE2D0FF","#1A9993FF"),
 UKBB = c(
   "#00441B", "#46A040", "#00AF99", "#FFC179", "#98D9E9", "#F6313E",
   "#8F1336", "#0081C9", "#001588", "#490C65", "#A65AC2", "#FF81AF",
   "#CDA2DB", "#FF5A00", "#FFCE00", "#FFA300", "#FF7900", "#FFA300" ),
 TF1 = c(
   "#F6313E", "#46A040", "#0081C9", "#A65AC2", "#FFA300",  "#FFFF32",
   "#89774A", "#FF6A80", "#999999", "#0DB2AA", "#001588", "#00441B",
   "#CDA2DB", "#98D9E9", "#FF9999", "#FFC966", "#A9FBA9" ),
 colx22 = c(
     '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
     '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
     '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
     '#000075', '#808080', '#4f34ff', '#f340F0'),
 paired = c(
     "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
     "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"),
 simpson = c("#FED439FF","#709AE1FF","#8A9197FF","#D2AF81FF","#FD7446FF",
             "#D5E4A2FF","#197EC0FF","#CC4C02","#46732EFF","#71D0F5FF",
             "#370335FF","#075149FF","#C80813FF","#91331FFF","#1A9993FF",
             "#FD8CC1FF"),
 tableau20 = c(
     "#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", "#98DF8A",
     "#D62728", "#FF9896", "#9467BD", "#C5B0D5", "#8C564B", "#C49C94",
     "#E377C2", "#F7B6D2", "#7F7F7F", "#C7C7C7", "#BCBD22", "#DBDB8D",
     "#17BECF", "#9EDAE5"),
 tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                     "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                     "#CDCC5D", "#6DCCDA"),
 colorblind10 = c("#006BA4", "#FF800E", "#ABABAB", "#595959",
                  "#5F9ED1", "#C85200", "#898989", "#A2C8EC",
                  "#FFBC79", "#CFCFCF"),
 trafficlight = c("#B10318", "#DBA13A", "#309343", "#D82526",
    "#FFC156", "#69B764", "#F26C64", "#FFDD71",
    "#9FCD99"),
 purplegray12 = c("#7B66D2", "#A699E8", "#DC5FBD", "#FFC0DA",
    "#5F5A41", "#B4B19B", "#995688", "#D898BA",
    "#AB6AD5", "#D098EE", "#8B7C6E", "#DBD4C5"),
 cyclic = c("#1F83B4", "#1696AC", "#18A188", "#29A03C", "#54A338",
    "#82A93F", "#ADB828", "#D8BD35", "#FFBD4C", "#FFB022",
    "#FF9C0E", "#FF810E", "#E75727", "#D23E4E", "#C94D8C",
    "#C04AA7", "#B446B3", "#9658B1", "#8061B4", "#6F63BB"),
 greenorange12 = c("#32A251", "#ACD98D", "#FF7F0F", "#FFB977",
                   "#3CB7CC", "#98D9E4", "#B85A0D", "#FFD94A",
                   "#39737C", "#86B4A9", "#82853B", "#CCC94D"),
alphabet = c(
  "#F0A0FF", "#0075DC", "#993F00", "#4C005C", "#191919", "#005C31",
  "#2BCE48", "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00",
  "#C20088", "#003380", "#FFA405", "#FFA8BB", "#426600", "#FF0010",
  "#5EF1F2", "#00998F", "#E0FF66", "#740AFF", "#990000", "#FFFF80",
  "#FFE100", "#FF5005"
),
alphabet2 = c(
  "#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656", "#1C8356",
  "#16FF32", "#F7E1A0", "#E2E2E2", "#1CBE4F", "#C4451C", "#DEA0FD",
  "#FE00FA", "#325A9B", "#FEAF16", "#F8A19F", "#90AD1C", "#F6222E",
  "#1CFFCE", "#2ED9FF", "#B10DA1", "#C075A6", "#FC1CBF", "#B00068",
  "#FBE426", "#FA0087"
),
glasbey = c(
  "#0000FF", "#FF0000", "#00FF00", "#000033", "#FF00B6", "#005300",
  "#FFD300", "#009FFF", "#9A4D42", "#00FFBE", "#783FC1", "#1F9698",
  "#FFACFD", "#B1CC71", "#F1085C", "#FE8F42", "#DD00FF", "#201A01",
  "#720055", "#766C95", "#02AD24", "#C8FF00", "#886C00", "#FFB79F",
  "#858567", "#A10300", "#14F9FF", "#00479E", "#DC5E93", "#93D4FF",
  "#004CFF", "#F2F318"
),
polychrome = c(
  "#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE",
  "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD",
  "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", "#85660D",
  "#B10DA1", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0",
  "#C075A6", "#782AB6", "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5",
  "#7ED7D1", "#1C7F93", "#D85FF7", "#683B79", "#66B0FF", "#3B00FB"
),
stepped = c(
  "#990F26", "#B33E52", "#CC7A88", "#E6B8BF", "#99600F", "#B3823E",
  "#CCAA7A", "#E6D2B8", "#54990F", "#78B33E", "#A3CC7A", "#CFE6B8",
  "#0F8299", "#3E9FB3", "#7ABECC", "#B8DEE6", "#3D0F99", "#653EB3",
  "#967ACC", "#C7B8E6", "#333333", "#666666", "#999999", "#CCCCCC" 
),
cellchat = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
    "#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", 
    "#222F75", "#1B9E77", "#B2DF8A", "#E3BE00", "#FB9A99", 
    "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261", 
    "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", "#B3DE69", 
    "#8DD3C7", "#999999"
    ),
cellchat_custom = c("#EC5757", "#377EB8", "#4DAF4A", "#984EA3", 
    "#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", 
    "#222F75", "#1B9E77", "#B2DF8A", "#E3BE00", "#FB9A99", 
    "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261", 
    "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", "#B3DE69", 
    "#8DD3C7", "#999999"
    ),
colx22_custom = c(
    '#EC5379', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
    '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
    '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
    '#000075', '#808080', '#4f34ff', '#f340F0'
    ),
cellpaper_custom = c(
    "#C1629D","#F36C65","#76B9CF","#C9A581","#EEB053","#94CFBF","#EAED70","#CD96BE","#CDDE8B","#E9B19C",
    "#ED4B43","#F39881","#AA4286","#3975E3","#77A879","#B589BC","#ECB888","#FC7440","#A032CB","#CD8AB1","#E9D6EB","#76A1EC","#A3C4A5","#EED785","#D4B595","#EFE0E7","#EEBED6","#F7D39B","#CD7560","#7FC28E","#F9CABE","#C76DA8","#FEF6AE","#6276B2","#E69D8E","#D9706F","#C5D9DF","#E16DB7","#908EBC","#AF88BB","#DEDBEE","#F07590","#DFB9D5","#B099B5","#5394C3","#C5B8B0","#917393","#97CADC","#DBEBF7","#84D7F7","#D4EEF6","#F5BCC9","#5EB4F3","#F7F4B7","#6F89D3","#50C2FB","#EA94CD","#D8B6F7","#EFC1EC"
    ),
cellpaper = c(
    "#D492BC","#F8ABA7","#A4D0DF","#DEC8B1","#F4CA8C","#C4E5DC","#F3F5AE","#E3C3DA","#E4EDBF","#F7E3DB",
    "#F48F89","#F9CABE","#C76DA8","#77A1EC","#A3C4A5","#FC7440","#A032CB","#CB9FD0","#EED785","#CBA57E","#CD7560","#7FC28E","#6276B2","#E16DB7","#908EBC","#AF88BB","#DEDBEE","#639EC8","#84D7F7","#917393"),

customecol2_light = c("#99D499FF", "#CBBEDDFF", "#FDCD9EFF", "#6089C0FF", "#F33599FF", "#B57162FF",
                    "#858585FF", "#49B192FF", "#E17F35FF", "#918DC2FF", "#D94955FF", "#69BDE5FF",
                    "#F2CA51FF", "#6D578FFF", "#D76891FF", "#599D33FF", "#E1DF62FF", "#BEDDFFFF",
                    "#F0C7EAFF", "#9A5A61FF", "#9065DBFF", "#696F98FF", "#EBA6A5FF", "#4D6E33FF",
                    "#F28133FF", "#97AD90FF", "#EB6385FF", "#D96DD2FF", "#FFE93BFF", "#B98BA9FF",
                    "#E19F6EFF", "#BDD16EFF", "#7685E5FF", "#DD69B4FF", "#76DE84FF", "#7EADD2FF",
                    "#338690FF", "#DBBD56FF", "#C4DFB7FF", "#B19D68FF", "#71455FFF", "#E99BB9FF",
                    "#8C6A4AFF", "#56D9E9FF", "#BE8190FF", "#B885CBFF", "#B55555FF", "#FFFFADFF",
                    "#65E1D5FF", "#CC7C45FF" ),
test1 = c(
    "#C4ABB1","#D29279","#8A569A","#F3A2BE","#D05E79","#00569A","#ECE05E","#ECC56E","#D5992B",
    "#70C6EB","#039DAE","#B2DDC4","#49986B","#A2C7C6","#2A86B7","#6364AF","#493F8B","#233671"),
pantone = c(
    "#5F4B8B", "#FF6F61", "#0F4C81", "#F5DF4D", "#9BB7D4", "#5ab481",
    "#be1456", "#E0B589", "#53B0AE", "#DD4124", "#3b1e8a", "#ec6887",
    "#855c4a", "#BB2649", "#028f6e", "#D5992B", "#88B04B",
    "#92A8D1", "#BCBABE", "#A89985", "#B565A7", "#D2386C",
    "#00A591", "#4B5338", "#EC5233", "#D29279", "#8A569A", 
    "#F3A2BE", "#D05E79", "#00569A", "#ECE05E", "#70C6EB", "#039DAE", 
    "#44398f"
),
group_default = c("#223D6C","#D20A13","#FFD121","#088247","#B037C4",
                  "#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D",
                  "#6E568C","#E0367A","#D8D155","#64495D","#7CC767"),
sampleid_default = c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#cea5c7',
                    '#53738c','#a5a9b0','#a78982','#696a6c','#92699e',
                    '#d69971','#df5734','#6c408e','#ac6894','#d4c2db',
                    '#537eb7','#83ab8e','#ece399','#405993','#cc7f73',
                    '#b95055','#d5bb72','#bc9a7f','#e0cfda','#d8a0c0',
                    '#e6b884','#b05545','#d69a55','#64a776','#cbdaa9',
                    '#efd2c9','#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
                    '#b85292','#6d6fa0','#8d689d','#c8c7e1','#d25774',
                    '#c49abc','#927c9a','#3674a2','#9f8d89','#72567a',
                    '#63a3b8','#c4daec','#61bada','#b7deea','#e29eaf',
                    '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
                    '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462'),
polychrome_cell = c("#b81316","#ea5c6f","#f7905a","#fe9e37","#e187cb","#fb948d",
                    "#e2b159","#ebed6f","#b2db87","#8acc72","#7ee7bb","#64cccf",
                    "#1f78b4","#a9dce6","#a48cbe","#e4b7d6","#6a3d9a")
)
##############################################
#' customized continious colors
#'
#' @export
continuous_palette <- list(
  spectral =  rev(RColorBrewer::brewer.pal(11, "Spectral")),
  viridis = viridis::viridis(11),
  inferno = viridis::inferno(11),
  magma = viridis::magma(11),
  cividis = viridis::cividis(11),
  rdwibl = c("#406AA8", "white", "#D91216"),
  orangepurple = c("#5A5B99", "white", "#F18704"),
  greenpurple = c("#095525","#57aa5e","#c1e6ba","#f4f1f4","#caafd4","#8a569c","#41004c"),
  jet = c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000" ),
  bluered12 = c("#2C69B0", "#B5C8E2", "#F02720", "#FFB6B0", "#AC613C","#E9C39B", "#6BA3D6", "#B5DFFD", "#AC8763", "#DDC9B4", "#BD0A36", "#F4737A"),
  GrandBudapest = c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236"),
  horizon = c('#000075', '#2E00FF', '#9408F7', '#C729D6', '#FA4AB5', '#FF6A95', '#FF8B74', '#FFAC53', '#FFCD32', '#FFFF60'),
  horizonExtra =c("#000436", "#021EA9", "#1632FB", "#6E34FC", "#C732D5", "#FD619D", "#FF9965", "#FFD32B", "#FFFC5A"),
  blueYellow = c( "#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3", "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D"),
  sambaNight = c( '#1873CC', '#1798E5', '#00BFFF', '#4AC596', '#00CC00', '#A2E700', '#FFFF00', '#FFD200', '#FFA500'),
  solarExtra = c( '#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D'),
  Moonrise1 = c("#F3DF6C", "#CEAB07", "#D5D5D3", "#24281A"),
  Royal1 = c("#899DA4", "#C93312", "#FAEFD1", "#DC863B"),
  fireworks = c("white", "#2488F0", "#7F3F98", "#E22929", "#FCB31A"),
  Moonrise2 = c("#798E87", "#C27D38", "#CCC591", "#29211F"),
  Cavalcanti = c("#D8B70A", "#02401B", "#A2A475", "#81A88D", "#972D15"),
  Royal2 = c("#9A8822", "#F5CDB4", "#F8AFA8", "#FDDDA0", "#74A089"),
  GrandBudapest2 = c("#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4"),
  Moonrise3 = c("#85D4E3", "#F4B5BD", "#9C964A", "#CDC08C", "#FAD77B"),
  Chevalier = c("#446455", "#FDD262", "#D3DDDC", "#C7B19C"),
  Zissou = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"),
  FantasticFox = c("#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20"),
  Darjeeling = c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"),
  Rushmore = c("#E1BD6D", "#EABE94", "#0B775E", "#35274A" ,"#F2300F"),
  BottleRocket = c("#A42820", "#5F5647", "#9B110E", "#3F5151", "#4E2A1E", "#550307", "#0C1707"),
  Darjeeling2 = c("#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#000000"),
  BottleRocket2 = c("#FAD510", "#CB2314", "#273046", "#354823", "#1E1E1E"),

  # Others
  lawhoops = c('#371377','#7700FF','#9E0142',"#CB2314",'#FF0080', '#DC494C',
               '#F88D51', "#FAD510","#FFFF5F",'#88CFA4','#238B45',"#02401B",
               "#0AD7D3", "dodgerblue", "#046C9A", "#273046", "#A2A475", "#354823", 'grey35', "#1E1E1E"),

  corona = c("#1f77b4","#d62728","#2ca02c","#ff7f0e","#9467bd","#8c564b","#e377c2","#7f7f7f",
           "#bcbd22","#17becf","#ad494a","#e7ba52","#8ca252","#756bb1","#636363","#aec7e8",
           "#ff9896","#98df8a","#ffbb78","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7","#dbdb8d",
           "#9edae5","#e7969c","#e7cb94","#c7e9c0","#de9ed6","#d9d9d9"),

  # Palettes from JDB
  citric = c('#03AB11', '#FFF301'),
  citrus = c('#15E602', '#FFFF33'),
  berry = c('#7700FF', '#FF0080'),
  forest = c('#E0E0E0', '#6BDA61', '#15E602'),
  white_mango = c('#FFFFFF', '#FF00FF', '#FF0000'),
  white_tango = c('#FFFFFF', '#FF00FF', '#9500FF'),
  white_grove = c('#FFFFFF', '#FFF301', '#BAE61D', '#01DF29'),
  white_jungle = c('#FFFFFF', '#FFF301', '#03AB11'),
  white_orange = c('#FFFFFF', '#FFF301', '#FF7700'),

  horizon = c('#000033', '#000075', '#0000B6', '#0000F8', '#2E00FF', '#6100FF', '#9408F7',
              '#C729D6', '#FA4AB5', '#FF6A95', '#FF8B74', '#FFAC53', '#FFCD32', '#FFEE11',
              '#FFFF60'),
  horizon_extra = c('#000033', '#000075', '#0000B6', '#0000F8', '#2E00FF', '#6100FF', '#9408F7', '#C729D6',
                    '#FA4AB5', '#FF6A95', '#FF8B74', '#FFAC53', '#FFCD32', '#FFEE11', '#FFFF60', '#FFFFFF'),
  wolfgang_basic = c('#FFFFD9', '#EDF8B1', '#C7E9B4', '#7FCDBB', '#41B6C4', '#1D91C0', '#225EA8', '#253494', '#081D58'),
  wolfgang_extra = c('#FFFFFF', '#FCFED3', '#E3F4B1', '#ABDEB6', '#60C1BF', '#2A9EC1', '#206AAD', '#243996', '#081D58'),
  solar_flare = c('#3361A5', '#2884E7', '#1BA7FF', '#76CEFF', '#FFFFFF', '#FFE060', '#FA8E24', '#DA2828', '#A31D1D'),
  solar_glare = c('#3361A5', '#2884E7', '#1BA7FF', '#76CEFF', '#FCFCFC', '#FFE060', '#FA8E24', '#DA2828', '#A31D1D'),
  solar_basic = c('#214B85', '#1873CC', '#1E90FF', '#00BFFF', '#ACD8E5', '#D2D2D2', '#FFD700', '#ED2C2C', '#A31D1D'),
  solar_extra = c('#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D'),
  solar_blues = c('#FCFCFC', '#C0E4FD', '#75CEFE', '#0CB9FF', '#1BA7FF', '#1E95FF', '#2884E7', '#3072C5', '#3361A5'),
  solar_rojos = c('#FCFCFC', '#FFEDB0', '#FFDF5F', '#FEC510', '#FA8E24', '#F14C2B', '#DA2828', '#BE2222', '#A31D1D'),
  samba_color = c('#1B85ED', '#1AA2F6', '#00BFFF', '#4AC596', '#00CC00', '#A7D400', '#FFD700', '#FFBE00', '#FFA500'),
  samba_night = c('#1873CC', '#1798E5', '#00BFFF', '#4AC596', '#00CC00', '#A2E700', '#FFFF00', '#FFD200', '#FFA500'),
  samba_light = c('#00D1FF', '#03AB11', '#FFF301'),

  dusk_dawn = c('#98ABC5', '#8D91AD', '#827896', '#775F80', '#6B476B', '#93575B', '#B8684A', '#DB7933', '#FF8C00'),
  dark_cyan = c('#000000', '#0E2824', '#014C44', '#0A5F4F', '#13725A', '#19997F', '#1EC0A6', '#19DFD2', '#00FFFF'),
  dark_blue = c('#000000', '#00171F', '#002F3F', '#00475F', '#005F7F', '#00779F', '#008FBF', '#00A7DF', '#00BFFF'),
  dark_citrus = c('#000000', '#22350F', '#3B680C', '#529111', '#6ABB15', '#74DD0F', '#7FFF00', '#ADF121', '#D1E131'),
  dark_violet = c('#000000', '#1E0A35', '#31016A', '#4B0181', '#660099', '#7800CA', '#8A00FF', '#C800FF', '#FE00FF'),
  ocean_green = c('#07519B', '#2975B4', '#5097C9', '#93C1DF', '#FCFCFC', '#CAEAC5', '#97D494', '#5BAB5A', '#006400'),
  ocean_earth = c('#0F3341', '#1563AA', '#0B99E6', '#3DCDFD', '#F7F7F7', '#B87350', '#872E1C', '#601622', '#401C2A'),
  ocean_brick = c('#0F3341', '#1563AA', '#0B99E6', '#3DCDFD', '#F7F7F7', '#EB9457', '#D1551F', '#B02F1B', '#8D1616'),
  algae_earth = c('#543005', '#985D12', '#CFA154', '#F0DEB1', '#F5F5F5', '#B5E2DC', '#5AB2A8', '#0E726A', '#003C30'),
  flame_flame = c('#000033', '#0000A5', '#1E00FB', '#6F00FD', '#C628D6', '#FE629D', '#FF9B64', '#FFD52C', '#FFFF5F'),
  flame_light = c('#000033', '#000E92', '#1300FF', '#8E0EEA', '#C628D6', '#E9699F', '#FF9B63', '#FFCE62', '#FFFF5F'),
  flame_polar = c('#C628D6', '#8E0EEA', '#1300FF', '#000E92', '#000033', '#7F494D', '#FF9B63', '#FFCE62', '#FFFF5F'),
  flame_volts = c('#000000', '#371377', '#5F00FF', '#9400FF', '#BE00FF', '#E000EB', '#FF00D8', '#FF0090', '#FF004B'),
  flame_watts = c('#FFFFFF', '#C190FF', '#5F00FF', '#9400FF', '#BE00FF', '#E000EB', '#FF00D8', '#FF0090', '#FF004B'),
  flame_artic = c('#000000', '#371377', '#5F00FF', '#BD00EC', '#FF00D8', '#C7ACEC', '#00FFFF', '#0AD7D3', '#0DB2AA'),
  flame_weird = c('#00FFFF', '#0AD7D3', '#0DB2AA', '#1C5551', '#000000', '#371377', '#5F00FF', '#BD00EC', '#FF00D8'),
  flame_blind = c('#0DB2AA', '#0AD7D3', '#00FFFF', '#B1FFFE', '#FFFFFF', '#FFA3EC', '#FF00D8', '#BD00EC', '#5F00FF'),
  flame_macaw = c('#000000', '#28410F', '#477C0E', '#64B114', '#9FCF23', '#C9E553', '#81F7D0', '#16DCD2', '#1AA58C'),
  flame_wings = c('#D1E131', '#85C51D', '#529111', '#2F4E0F', '#000000', '#0F4338', '#107E6A', '#1BBBA7', '#00FFFF'),

  calma_azules = c('#031C25', '#093B4D', '#1C5F77', '#3685A2', '#56A6C3', '#86C2D8', '#B6DDEB', '#F2FBFE'),
  calma_musgos = c('#212503', '#444D09', '#6B771C', '#93A236', '#B4C356', '#CDD886', '#E4EBB6', '#FCFEF2'),
  calma_bosque = c('#032506', '#094D0E', '#1C7722', '#36A23D', '#56C35D', '#86D88B', '#B6EBBA', '#F2FEF3'),
  calma_marino = c('#032515', '#094D2D', '#1C774D', '#36A26F', '#56C390', '#86D8B2', '#B6EBD2', '#F2FEF8'),
  calma_morado = c('#030925', '#09154D', '#1C2B77', '#3648A2', '#5668C3', '#8694D8', '#B6BFEB', '#F2F4FE'),
  calma_manudo = c('#290303', '#590707', '#8C1616', '#BE2A2A', '#DF4A4A', '#ED8080', '#F7B4B4', '#FFEEEE'),
  china_theory = c('#120324', '#420A4A', '#721D57', '#9B3850', '#BC6B58', '#D3B687', '#E6E8B7', '#F8FDF2'),
  china_ranges = c('#031424', '#1F0A4A', '#721D64', '#9B3838', '#BCAB58', '#A0D387', '#B7E8CF', '#F2F9FD'),
  china_weirdo = c('#04032E', '#2E0267', '#890BA3', '#DE15AF', '#FF347E', '#FF7772', '#FFCFAB', '#FFFBEA'),
  china_basics = c('#25032E', '#670253', '#A30B48', '#DE1515', '#FF8534', '#FFE272', '#EEFFAB', '#F2FFEA'),
  china_sunset = c('#031124', '#0C0A4A', '#451D72', '#91389B', '#BC589B', '#D38799', '#E8C1B7', '#FDF9F2'),
  china_dragon = c('#03032A', '#2B065C', '#801491', '#C52696', '#E74671', '#F2917D', '#FADEB3', '#FEFFED'),
  china_novice = c('#2A0E03', '#5C4406', '#7E9114', '#68C526', '#46E748', '#7DF2B2', '#B3FAF2', '#EDF9FF'),

  brewer_fire = c('#FFFFE5', '#FFF7BC', '#FEE391', '#FEC44F', '#FE9929', '#EC7014', '#CC4C02', '#993404', '#662506'),
  brewer_heat = c('#FFF7EC', '#FEE8C8', '#FDD49E', '#FDBB84', '#FC8D59', '#EF6548', '#D7301F', '#B30000', '#7F0000'),
  brewer_orange = c('#FFF5EB', '#FEE6CE', '#FDD0A2', '#FDAE6B', '#FD8D3C', '#F16913', '#D94801', '#A63603', '#7F2704'),
  brewer_red = c('#FFF5F0', '#FEE0D2', '#FCBBA1', '#FC9272', '#FB6A4A', '#EF3B2C', '#CB181D', '#A50F15', '#67000D'),
  brewer_green = c('#F7FCF5', '#E5F5E0', '#C7E9C0', '#A1D99B', '#74C476', '#41AB5D', '#238B45', '#006D2C', '#00441B'),
  brewer_blue = c('#F7FBFF', '#DEEBF7', '#C6DBEF', '#9ECAE1', '#6BAED6', '#4292C6', '#2171B5', '#08519C', '#08306B'),
  brewer_purple = c('#FCFBFD', '#EFEDF5', '#DADAEB', '#BCBDDC', '#9E9AC8', '#807DBA', '#6A51A3', '#54278F', '#3F007D'),
  brewer_violet = c('#FFF7F3', '#FDE0DD', '#FCC5C0', '#FA9FB5', '#F768A1', '#DD3497', '#AE017E', '#7A0177', '#49006A'),
  brewer_jamaica = c('#006837', '#2DA154', '#86CB66', '#CCE982', '#FFFFBF', '#FDD380', '#F88D51', '#DE3F2E', '#A50026'),
  brewer_marine = c('#F7FCF0', '#E0F3DB', '#CCEBC5', '#A8DDB5', '#7BCCC4', '#4EB3D3', '#2B8CBE', '#0868AC', '#084081'),
  brewer_spectra = c('#5E4FA2', '#3F96B7', '#88CFA4', '#D7EF9B', '#FFFFBF', '#FDD380', '#F88D51', '#DC494C', '#9E0142'),
  brewer_celsius = c('#313695', '#5083BB', '#8FC3DD', '#D2ECF4', '#FFFFBF', '#FDD384', '#F88D51', '#DE3F2E', '#A50026'),
  brewer_yes = c('#053061', '#2971B1', '#6AACD0', '#C1DDEB', '#F7F7F7', '#FACDB5', '#E58267', '#BB2933', '#67001F'),

  forest_yellow = c('#215e33', '#306835', '#3f7136', '#4e7b38', '#5d843a', '#6c8e3b', '#7b983d', '#89a13f', '#98ab41', '#a7b542', '#b6be44', '#c5c846', '#d4d147', '#e3db49'),
  forest_citric = c('#0b3310', '#1b4210', '#2b5111', '#3b6111', '#4c7012', '#5c7f12', '#6c8e12', '#7c9e13', '#8cad13', '#9cbc13', '#adcb14', '#bddb14', '#cdea15', '#ddf915'),
  citric_yellow = c('#21be45', '#30c246', '#40c546', '#4fc947', '#5fcc48', '#6ed048', '#7dd349', '#8dd74a', '#9cda4b', '#abde4b', '#bbe14c', '#cae54d', '#dae84d', '#e9ec4e'),
  ocean_citrus  = c('#3683ba', '#418bb0', '#4c93a7', '#569b9d', '#61a393', '#6cab8a', '#77b380', '#81ba76', '#8cc26c', '#97ca63', '#a2d259', '#acda4f', '#b7e246', '#c2ea3c'),
  ocean_pink = c('#3b5e84', '#4a5e84', '#595e84', '#685e84', '#765e84', '#855e84', '#945e84', '#a35f85', '#b25f85', '#c15f85', '#cf5f85', '#de5f85', '#ed5f85', '#fc5f85'),
  ocean_red = c('#3b82ae', '#487aa1', '#547294', '#616987', '#6d617a', '#7a596d', '#865160', '#934853', '#9f4046', '#ac3839', '#b8302c', '#c5271f', '#d11f12', '#de1705'),
  ocean_aqua = c('#2668aa', '#2d6eaa', '#3574aa', '#3c7aaa', '#4380a9', '#4b86a9', '#528ca9', '#5993a9', '#6099a9', '#689fa9', '#6fa5a8', '#76aba8', '#7eb1a8', '#85b7a8'),
  ocean_teal = c('#0a0a66', '#15176c', '#202573', '#2b3279', '#364080', '#414d86', '#4c5a8d', '#576893', '#62759a', '#6d82a0', '#7890a7', '#839dad', '#8eabb4', '#99b8ba'),
  cyan_brick = c('#6cd0c2', '#70c3b3', '#74b6a5', '#78a996', '#7b9c88', '#7f8f79', '#83826a', '#87755c', '#8b684d', '#8f5b3e', '#924e30', '#964121', '#9a3413', '#9e2704'),
  aqua_brick = c('#019bcf', '#0d94c2', '#1a8cb4', '#2685a7', '#337d99', '#3f768c', '#4c6e7e', '#586771', '#655f63', '#715856', '#7e5048', '#8a493b', '#97412d', '#a33a20'),
  aqua_tan  = c('#62adbb', '#6cb0b5', '#75b2af', '#7fb5a8', '#88b7a2', '#92ba9c', '#9bbd96', '#a5bf8f', '#aec289', '#b8c583', '#c1c77d', '#cbca76', '#d4cc70', '#decf6a'),
  cyan_tan  = c('#4fe8c2', '#5be3bb', '#67deb5', '#73d9ae', '#7fd3a7', '#8bcea1', '#97c99a', '#a4c493', '#b0bf8c', '#bcba86', '#c8b47f', '#d4af78', '#e0aa72', '#eca56b'),
  teal_orange = c('#0cb499', '#1dae92', '#2ea98b', '#3fa384', '#4f9e7d', '#609876', '#71936f', '#828d69', '#938862', '#a4825b', '#b47d54', '#c5774d', '#d67246', '#e76c3f'),
  teal_violet = c('#4b8b84', '#508184', '#567784', '#5b6d84', '#616384', '#665984', '#6b4f84', '#714683', '#763c83', '#7b3283', '#812883', '#861e83', '#8c1483', '#910a83'),
  blue_cyan = c('#4111f2', '#4923f0', '#5135ed', '#5946eb', '#6258e8', '#6a6ae6', '#727ce3', '#7a8de1', '#829fde', '#8ab1dc', '#93c3d9', '#9bd4d7', '#a3e6d4', '#abf8d2'),
  purple_pink = c('#6848d1', '#7345ca', '#7f42c3', '#8a3fbc', '#953db5', '#a13aae', '#ac37a7', '#b734a1', '#c2319a', '#ce2e93', '#d92c8c', '#e42985', '#f0267e', '#fb2377'),
  purple_baby = c('#511293', '#5a239b', '#6234a2', '#6b45aa', '#7356b2', '#7c67b9', '#8478c1', '#8d88c9', '#9599d1', '#9eaad8', '#a6bbe0', '#afcce8', '#b7ddef', '#c0eef7'),
  cyan_green = c('#29ddea', '#31dcd8', '#39dcc7', '#41dbb5', '#4adba3', '#52da92', '#5ad980', '#62d96e', '#6ad85c', '#72d74b', '#7bd739', '#83d627', '#8bd616', '#93d504'),
  cyan_pink = c('#606bef', '#6c6ae6', '#7869dd', '#8468d4', '#9168cb', '#9d67c2', '#a966b9', '#b565b1', '#c164a8', '#cd639f', '#da6396', '#e6628d', '#f26184', '#fe607b'),
  cyan_violet = c('#29e9ae', '#36dbae', '#43cdae', '#50beaf', '#5eb0af', '#6ba2af', '#7894af', '#8585b0', '#9277b0', '#9f69b0', '#ad5bb0', '#ba4cb1', '#c73eb1', '#d430b1'),
  cyan_purple = c( '#4adbf1', '#51cff0', '#59c3ef', '#60b7ee', '#68abed', '#6f9fec', '#7793eb', '#7e88eb', '#867cea', '#8d70e9')
)
## 提取颜色
get_colors <- function(object_meta, groupby, palette = "customecol2", start_from = NULL){

                groupby_col = paste0(groupby,"_col")
                if ((! groupby_col %in% colnames(object_meta)) & groupby %in% colnames(object_meta)){
                    print(paste0(groupby,"颜色信息不存在，写入",groupby,"_col"))
                    # 适配多种颜色选择
                    group_colors = list()
                    for ( x in  unlist(strsplit(palette, ",", perl =T))){
                        m = unlist( strsplit( x, ":", perl =T))
                        if(length(m) == 1 ){
                            group_colors[[groupby]] = m
                        }else if(length(m) == 2){
                            if (m[1] == groupby){
                                group_colors[[m[1]]] = m[2]
                                print(paste0("指定",m[1],"列色板为",m[2],"，不采用原注释颜色（如有）。"))
                                object_meta = object_meta[ ,!grepl(paste0("^",m[1],"_col$" ), colnames(object_meta))]
                            }
                        }else{
                            stop(paste0(palette, " 填写有误，应为：groupby1:palette1,groupby2:palette2 格式。"))
                        }
                    }
                    #当start_from未被定义，且包含关键词时，自动从第11位截取颜色
                    if (is.null(start_from)){
                        #if(grepl("^group", groupby) || groupby == "sampleid") {
                            #print(paste0('start_from未被赋值且分组变量',groupby,'为sampleid或以group开头，将自动从11位开始取色'))
                            #start_from = 11
                        #}else{
                            start_from = 1
                        #}
                    }
                    if( groupby %in% names(group_colors) ){
                        palette_use = group_colors[[groupby]]
                    }else{
                        print(paste0("rds中无 ",groupby_col," 颜色列且色板未指定，采用默认customecol2色板!"))
                        palette_use = "customecol2"
                    }
                    if ( !is.null(discrete_palette[[palette_use]]) ){
                        # 从第{start_from}位色板开始取色
                        palettes <- discrete_palette[[palette_use]][start_from:length(discrete_palette[[palette_use]])]
                        if( length(palettes) < length(unique(object_meta[,groupby]))){
                            print(paste0(groupby,"的色板不够用了，自动生成颜色，如需指定颜色，请提供color_file文件!"))
                            palettes = grDevices::colorRampPalette(palettes)(length(unique(object_meta[,groupby])))
                        }
                    }else{
                        stop(paste0("rds中无 ",groupby_col," 颜色列且色板填写有误，请核查!"))
                    }
                    if(class(object_meta[,groupby]) == 'factor'){
                        all_types = levels(object_meta[,groupby])
                    } else {
                        all_types = sort(unique(object_meta[,groupby]))
                    }
                    col_tmp = palettes[1:length(all_types)]
                    i = 1
                    object_meta[,groupby_col]=""
                    for ( celltype in all_types ){
                      object_meta[which(object_meta[,groupby] == celltype),groupby_col] = as.character(col_tmp[i])
                      i = i + 1
                    }
                }
                print(paste0("根据 ",groupby_col," 颜色列取色。"))
                if ( length(unique(object_meta[,groupby])) != length(unique(object_meta[,groupby_col]))){
                    stop(paste0(groupby,"颜色数量",length(unique(object_meta[,groupby_col])),"与细胞群数量",length(unique(object_meta[,groupby])),"不对应！请核查rds！"))
                }
                if ( class(object_meta[,groupby])=="factor" ){
                    nlevel_list = levels(object_meta[,groupby])
                }else{
                    nlevel_list = sort(as.character(unique(object_meta[,groupby])))
                }
                tmp_df <- unique(object_meta[,c(groupby,groupby_col)])
                new_celltype_pal <- as.vector(tmp_df[,groupby_col])
                names(new_celltype_pal) <-  as.vector(tmp_df[,groupby])
                new_celltype_pal = new_celltype_pal[nlevel_list]
                user_color_pal = as.vector(new_celltype_pal)
                print("对应颜色选择为：")
                print(new_celltype_pal)
                result = list(user_color_pal = user_color_pal, new_celltype_pal = new_celltype_pal, object_meta = object_meta)
                return(result)
            }
## 注释颜色
color_anno <- function(object_meta, color_file){
        groupby_all = colnames(color_file)[!grepl("_col$", colnames(color_file))]
        print(paste0("传入文件中共检测到",length(groupby_all),"列元素需要注释颜色"))
        for (groupby in groupby_all){
            print(paste0("对",groupby,"列元素注释颜色"))
            #anno_groupby = colnames(color_file)[1]
            groupby_col = paste0(groupby,"_col")
            if( !groupby %in% colnames(object_meta)){
                stop(paste0("输入文件",groupby,"列在rds中不存在，请检查输入文件！"))
            }
            color_match = as.vector(color_file[,groupby_col])
            names(color_match) = as.vector(color_file[,groupby])
            color_match = color_match[ !is.na( match(names(color_match), object_meta[,groupby])) ]
            if( length(setdiff( names(color_match), unique(as.vector(object_meta[,groupby])) )) != 0 ){
                stop(paste0("输入文件 ",groupby," 列元素与rds中该列元素长度不等，请检查输入文件！"))
            }
            object_meta[, groupby_col] = color_match[as.vector(object_meta[,groupby])]
        }
        return(object_meta)
        }

get_color_use <- function(data_ob, groupby, palette = "customecol2", start_from = NULL, color_file = NULL, use_color_anno = T){
        print(paste0("开始对",groupby,"取色,色板为：",palette,"；是否采用seurat对象中的颜色注释：", use_color_anno))
        if ( ! as.logical(use_color_anno )){
            print(paste0("删除原",groupby,"_col列颜色注释重新取色。"))
            data_ob@meta.data = data_ob@meta.data[ ,!grepl(paste0("^",groupby,"_col$" ), colnames(data_ob@meta.data))]
        }
        if ( !is.null(color_file)){
            print(paste0("根据颜色文件 ",color_file," 注释颜色。"))
            color_file = read.delim(color_file, sep="\t", header = T)
            meta_anno = color_anno(data_ob@meta.data, color_file)
        } else {
            meta_anno = data_ob@meta.data
        }
        color_use = get_colors(meta_anno, groupby, palette, start_from)
        return(color_use)
        }



