
String filename = "fmm_data_ellipse_elem_50.conf";
String elem_val_in_file = "";
int imgWidth = 600;
String output_file = "output.png";
final String[] types = {"POINTS","ELEM"};

String type = "ELEM";
boolean values = elem_val_in_file.length() > 0;


String[] input;
Pot_Elem[] elem;
float[] min_max_node;

Boolean isValidType(String t)
{
  Boolean isValid = false;
  for(int i = 0; i<types.length ; i++)
  {
    if (t.equals(types[i]))
      isValid = true;
  }
  return isValid;
}

void drawPoints(String[] conf, int imgWidth, String out_file_name)
{
  //set colors for source and target points
  color source_col = #FF0000; //red
  color tgt_col = #0000FF; //blue
  color tgt_src_col = #FF00FF; //violet
  Pot_Point[] points = new Pot_Point[conf.length-1];
  println("number of points " + str(points.length));
  
  float maxX = 0;
  float maxY = 0;
  float minX = 10000000;
  float minY = 10000000;
  
  //read the points and get max X and max Y
  for(int i = 1; i<conf.length; i++)
  {
    String line = conf[i];
    String[] tokens = split(line,' ');
    if(tokens.length < 4)
    {
      println("there was an error parsing the data");
      exit();
    }
    float posX = float(tokens[0]);
    float posY = float(tokens[1]);
    if(maxX < posX) maxX = posX;
    if(maxY < posY) maxY = posY;
    if(minX > posX) minX = posX;
    if(minY > posY) minY = posY;
    int source = int(tokens[2]);
    int target = int(tokens[3]);
    int val = 2*target + source;
    color point_col;
    
    switch(val) {
      case 1:
        point_col = source_col;
        break;
      case 2:
        point_col = tgt_col;
        break;
      case 3:
        point_col = tgt_src_col;
        break;
      default:
        point_col = color(#FFFFFF);
        println("There is a point that's neither source or target");
        exit();
    }
    points[i-1] = new Pot_Point(posX,posY,point_col); //<>//
  }
  
  float xrange = maxX-minX;
  float yrange = maxY-minY;
  
  // create Image with given size
  int imgHeight = int(yrange/xrange*imgWidth);
  PImage img = createImage(imgWidth, imgHeight, RGB);
  img.loadPixels();
  for(int i = 0; i<img.pixels.length; i++)
    img.pixels[i] = color(#FFFFFF);
    
  // set pixels
  for(int i = 0; i< points.length; i++)
  {
    int x = int((points[i].x-minX) * (imgWidth-1)/xrange);
    int y = imgHeight-1 - int((points[i].y-minY) * (imgHeight-1)/yrange);
    img.pixels[imgWidth*y + x] = points[i].col;
  }
  
  img.updatePixels();
  img.save(out_file_name);
}

Pot_Elem[] readElements(String[] conf)
{
   //set colors for source element and target points
  color source_col = #FF0000; //red
  color tgt_col = #0000FF; //blue
  
  // get number of elements
  int num_el = int(conf[1]);
  // save nodes in array
  float[][] nodes = new float[num_el][];
  Pot_Elem[] elems = new Pot_Elem[num_el];
  for(int i = 0; i<num_el;i++)
  {
    String line = conf[2+i];
    String[] tokens = split(line,' ');
    float x = float(tokens[1]);
    float y = float(tokens[2]);
    
    nodes[i] = new float[2];
    nodes[i][0] = x;
    nodes[i][1] = y;
  }
  
  // read element configuration
  for(int i = 2+num_el; i<2+2*num_el; i++)
  {
    String line = conf[i];
    String[] tokens = split(line,' ');
    int startNode = int(tokens[0]);
    int endNode = int(tokens[1]);
    int src = int(tokens[2]);
    int tgt = int(tokens[3]);
    elems[i-2-num_el] = new Pot_Elem(nodes[startNode][0],
                                     nodes[startNode][1],
                                     nodes[endNode][0],
                                     nodes[endNode][1],
                                     boolean(src),
                                     boolean(tgt));
   if(boolean(src)) elems[i-2-num_el].src_col = source_col;
   if(boolean(tgt)) elems[i-2-num_el].tgt_col = tgt_col;
   }
   
   return elems;
}

//0:maxX, 1:maxY, 2:minX, 3:minY
float[] min_max_nodes(Pot_Elem[] elems)
{
  float maxX = 0;
  float maxY = 0;
  float minX = 10000000;
  float minY = 10000000;
  
  //read the points and get max X and max Y
  for(int i = 1; i<elems.length; i++)
  {
    float posX1 = max(elems[i].startX,elems[i].endX);
    float posY1 = max(elems[i].startY,elems[i].endY);
    float posX2 = min(elems[i].startX,elems[i].endX);
    float posY2 = min(elems[i].startY,elems[i].endY);
    if(maxX < posX1) maxX = posX1;
    if(maxY < posY1) maxY = posY1;
    if(minX > posX2) minX = posX2;
    if(minY > posY2) minY = posY2;
  }
  float[] res = {maxX,maxY,minX,minY};
  return res;
}

void draw_Elements(float[] min_max, Pot_Elem[] elem, String out_file, float min_val, float max_val)
{  
  float xrange = min_max[0]-min_max[2];
  float yrange = min_max[1]-min_max[3];
  float minX = min_max[2];
  float minY = min_max[3];
  
  int h = int(imgWidth*yrange/xrange);
  int barWidth = imgWidth/2;
  int barDist = 0;
  int barHeight = 0;
  if(values && (max_val - min_val > 0))
  {
    // add a bar
    barDist = 30;
    barHeight = 10;
  }
  PGraphics pg = createGraphics(imgWidth,h+barDist+barHeight);
  pg.beginDraw();
  pg.colorMode(RGB,255);
  pg.background(255);
  
  if(values && (max_val - min_val > 0))
  {
    //draw bar
    pg.colorMode(HSB,360);
    
    for(int i =imgWidth/4; i<(3*imgWidth)/4;i++)
    {
      for(int j = h+barDist; j<h+barDist+barHeight; j++)
      {
        pg.stroke(240*(2*i-imgWidth/2)/imgWidth,360,360);
        pg.point(i,j);
      }
    }
    
    pg.colorMode(RGB,255);
  }
  
  for(int i = 0; i< elem.length; i++)
  {
    int start_x = int((elem[i].startX-minX) * (imgWidth-1)/xrange);
    int start_y = h-1 - int((elem[i].startY-minY) * (h-1)/yrange);
    int end_x = int((elem[i].endX-minX) * (imgWidth-1)/xrange);
    int end_y = h-1 - int((elem[i].endY-minY) * (h-1)/yrange);
    int center_x = int((elem[i].centerX-minX) * (imgWidth-1)/xrange);
    int center_y = h-1 - int((elem[i].centerY-minY) * (h-1)/yrange);
    if(!elem[i].hasValue)
    {
      if(elem[i].source) {
        pg.stroke(elem[i].src_col);
        pg.line(start_x,start_y,end_x,end_y);
      }
      if(elem[i].target) {
        pg.stroke(elem[i].tgt_col);
        pg.strokeWeight(3);
        pg.point(center_x,center_y);
        pg.strokeWeight(1);
      }
    } 
    else
    {
      int hue = 0;
      if(min_val < max_val)
        hue = int((elem[i].value - min_val)/(max_val-min_val)*240);
        pg.colorMode(HSB,360);
        pg.stroke(hue,360,360);
        pg.strokeWeight(2);
        pg.line(start_x,start_y,end_x,end_y);
        pg.strokeWeight(1);
        pg.colorMode(RGB,255);
    }
  }
  pg.endDraw();
  pg.save(out_file);
}

// 0:minValue, 1:maxValue
float[] setValues(Pot_Elem[] elem, String filename)
{
  float minVal, maxVal;
  float [] res = new float[2];
  String[] vals = loadStrings(filename);
  if(vals.length != elem.length)
  {
    println("number of elements is not number of values");
    exit();
  }
  minVal = float(vals[0]);
  maxVal = minVal;
  for(int i = 0; i<vals.length; i++)
  {
    String line = vals[i];
    float val = float(line);
    if(val > maxVal) maxVal = val;
    if(val < minVal) minVal = val;
    elem[i].value = val;
    elem[i].hasValue = true;
  }
  res[0] = minVal;
  res[1] = maxVal;
  println("minimal value = " + str(minVal));
  println("maximal value = " + str(maxVal));
  return res;
}


void setup() 
{
  size(600, 360);
  noLoop();
}

void draw()
{  
  background(255);
  input = loadStrings(filename);
  if(input == null)
  {
    println("Invalid input file, exiting now");
    exit();
  }
  if (!isValidType(type)) {
    println("Invalid type");
    exit();
  }
  //set window appropriate
  if(type.equals("ELEM"))
  {
    elem = readElements(input);
    min_max_node = min_max_nodes(elem);
  } 
  else
  {
    background(100);
  }
  if(type.equals("POINTS"))
  {
    drawPoints(input, imgWidth, output_file);
    exit();
  } else if (type.equals("ELEM"))
  {
    float min_val = 1;
    float max_val = 0;
    if(values)
    {
      float[] vals = setValues(elem,elem_val_in_file);
      min_val = vals[0];
      max_val = vals[1];
    }
    draw_Elements(min_max_node,elem, output_file, min_val, max_val);
    exit();
  }
}