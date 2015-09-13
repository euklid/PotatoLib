public class Pot_Point
{
  public Pot_Point()
  {
    x = 0;
    y = 0;
    col = #FFFFFF;  
  };
  
  public Pot_Point(float x, float y, color col)
  {
    this.x = x;
    this.y = y;
    this.col = col;
  }
  
  public float x,y;
  public color col;
}