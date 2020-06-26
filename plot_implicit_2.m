function  plot_implicit_2(c)
% c2=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
% c=c.*c2;


%[x,y]=meshgrid(-3:0.01:3);

syms y x

func=c(1)+c(2)*power(x,1)+c(3)*y+c(4)*power(x,2)+c(5)*x*y+c(6)*power(y,2);


ezplot(func,[-1,1]);

% package com.cammaway.garden_test;
% 
% import android.support.v7.app.AppCompatActivity;
% import android.os.Bundle;
% import android.view.KeyEvent;
% import android.view.View;
% import android.widget.AdapterView;
% import android.widget.ArrayAdapter;
% import android.widget.EditText;
% import android.widget.ListView;
% import android.widget.Spinner;
% import android.widget.TextView;
% import android.widget.Toast;
% 
% import java.util.ArrayList;
% import java.util.Arrays;
% 
% class ScoreSet{
%     int vegsInList;
%     String Vegs[]=new String[8];
%     String ScoreName[]={"HPscore","NITROGENscore"};
%     int scores[]=new int[8];
%     public void ScoreSet(String _Vegs)
%     {
%         vegsInList=0;
%     }
%     
%     public void addVeg(String _Vegs)
%     {
%         Vegs[vegsInList]=_Vegs;
%         vegsInList++;
%         assignScores();
%         
%     }
%     
%     public void assignScores()
%     {
%         for (int i=0;i<Vegs.length;i++)//go over all vegs
%         {
%             for (int j=0;j<ScoreName.length;j++)///assign scores
%             {
%                 scores[j]=5;
%                 
%             }
%             
%             
%         }
%     }
%     
%     
% }
% public class MainActivity extends AppCompatActivity {
%     Spinner spinner1;
%     ScoreSet Scores=new ScoreSet();
%     ArrayAdapter<String> adapterBusinessType;
%     String vegType[] = { "Carrot", "ChickPea", "Paw paw", "Beets",
%         "Sweet Potato"};
%     @Override
%     public void onCreate(Bundle savedInstanceState) {
%         super.onCreate(savedInstanceState);
%         // Inflate your View setContentView(R.layout.main);
%         // Get references to UI widgets
%         // Inflate your View
%         Button myButton = new Button(this);
%         RelativeLayout myLayout = new RelativeLayout(this);
%         myLayout.addView(myButton);
%         setContentView(myLayout);
%         // setContentView(R.layout.activity_main);
%         final EditText myEditText = new EditText(this);// findViewById(R.id.myEditText);
%         
%         final ListView myListView=(ListView) findViewById(R.id.myListView);
%         final TextView[] myTextViews=new TextView[2];
%         
%         for (int i=0;i<Scores.ScoreName.length;i++) {
%             
%             
%             myTextViews[i]=(TextView)findViewById(R.id.PHscore);
%         }
%         /* spinner1 = (Spinner) findViewById(R.id.vegList);
%          
%          adapterBusinessType= new ArrayAdapter(this, android.R.layout.simple_spinner_item,vegType);
%          spinner1.setAdapter(adapterBusinessType);
%          
%          spinner1.setOnItemSelectedListener(new AdapterView.OnItemSelectedListener() {
%          
%          @Override
%          public void onItemSelected(AdapterView<?> adapter, View v,
%          int position, long id) {
%          // On selecting a spinner item
%          String item = adapter.getItemAtPosition(position).toString();
%          
%          // Showing selected spinner item
%          Toast.makeText(getApplicationContext(),
%          "Selected Country : " + item, Toast.LENGTH_LONG).show();
%          }
%          
%          @Override
%          public void onNothingSelected(AdapterView<?> arg0) {
%          // TODO Auto-generated method stub
%          
%          }
%          });
%          */
%         
%         // Create the Array List of to do items
%         final ArrayList<String> todoItems = new ArrayList<String>();
%         // Create the Array Adapter to bind the array to the List View final ArrayAdapter<String> aa;
%         final ArrayAdapter<String> aa;
%         aa = new ArrayAdapter<String>(this, android.R.layout.simple_list_item_activated_1, todoItems);
%         // Bind the Array Adapter to the List View
%         myListView.setAdapter(aa);
%         
%         myEditText.setOnKeyListener(new View.OnKeyListener() {
%             public boolean onKey(View v, int keyCode, KeyEvent event) {
%                 if (event.getAction() == KeyEvent.ACTION_DOWN)
%                     if ((keyCode == KeyEvent.KEYCODE_DPAD_CENTER) ||
%                         (keyCode == KeyEvent.KEYCODE_ENTER)) {
%                         String new_veg=myEditText.getText().toString();
%                         int found=0;
%                         for (int i=0;i<vegType.length;i++)
%                         {
%                             if (vegType[i].equals(new_veg))//match found
%                             {
%                                 if (Arrays.asList(Scores.Vegs).contains(new_veg))
%                                 {
%                                     Toast.makeText(getApplicationContext(),
%                                                    "The entered vegtable already exists ", Toast.LENGTH_LONG).show();
%                                 }
%                                 else {
%                                     
%                                     
%                                     Scores.addVeg(new_veg);
%                                     ProcessVegInfo(new_veg);
%                                     todoItems.add(0, myEditText.getText().toString());
%                                     aa.notifyDataSetChanged();
%                                     myEditText.setText("");
%                                     for (int j=0;j<=Scores.ScoreName.length;j++) {
%                                         PHscore_text.setText("Score on " + Scores.ScoreName[j] + ":" + String.valueOf(Scores.scores[0]));
%                                     }
%                                     found = 1;
%                                     break;
%                                 }
%                                 
%                             }
%                             
%                         }
%                         if (found==0) {
%                             Toast.makeText(getApplicationContext(),
%                                            "The entered vegtable doesn't exist ", Toast.LENGTH_LONG).show();
%                         }
%                         return true;
%                     }
%                 return false;
%             } });
%     }
%     
%     private void ProcessVegInfo(String Veg)
%     {
%         
%         
%     }
%     
% }
