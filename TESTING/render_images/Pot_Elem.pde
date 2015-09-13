public class Pot_Elem
{
  public Pot_Elem(){
    this.hasValue = false;
  }
  public Pot_Elem(float startX, float startY, float endX, float endY, boolean src, boolean target)
  {
    this.startX = startX;
    this.startY = startY;
    this.endX = endX;
    this.endY = endY;
    this.hasValue = false;
    this.value = 0;
    this.centerX = (startX + endX)/2;
    this.centerY = (startY+endY)/2;
    this.source = src;
    this.target = target;
  }
  
  public Pot_Elem(float startX, float startY, float endX, float endY, boolean src, boolean target, float val)
  {
    this(startX,startY,endX,endY, src, target);
    this.hasValue = true;
    this.value = val;
  }
  
  public float startX,startY,endX,endY,centerX,centerY;
  public float value;
  public boolean hasValue;
  public boolean source, target;
  public color src_col, tgt_col;
}