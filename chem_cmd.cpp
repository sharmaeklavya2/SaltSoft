#include<iostream>
#include<string>
#include<conio.h>	//for _get_restrict_che()

void disp_error(const char* str)
{
#ifdef _WINDOWS_
	MessageBox(NULL,str,"Chem Analyze error",MB_ICONERROR|MB_OK);
#else
	std::cout<<"Error:\t"<<str<<std::endl;
#endif
}

void disp_error(const std::string& str)
{
#ifdef _WINDOWS_
	MessageBox(NULL,str.c_str(),"Chem Analyze error",MB_ICONERROR|MB_OK);
#else
	std::cout<<"Error:\t"<<str<<std::endl;
#endif
}

#include"chem.h"

char _get_restrict_che(const char* a)
{
	char ch;unsigned i;
	while(true)
	{
		ch=_getch();
		for(i=0;a[i]!='\0';i++)if(ch==a[i])
		{_putch(ch);return ch;}
	}
}

void disp_msg(const char* str)
{
#ifdef _WINDOWS_
	MessageBox(NULL,str,"Chem Analyze msg",MB_ICONINFORMATION|MB_OK);
#else
	std::cout<<str<<std::endl;
#endif
}

void disp_msg(const std::string& str)
{
#ifdef _WINDOWS_
	MessageBox(NULL,str.c_str(),"Chem Analyze msg",MB_ICONINFORMATION|MB_OK);
#else
	std::cout<<str<<std::endl;
#endif
}

using namespace std;

void reaction()
{
	std::string str;
	std::set<std::string> obs_set;
	std::set<std::string>::iterator it;
	tube_t tube;
	bool heat;
	do
	{
		cout<<"\nEnter reactants:\t";
		if(cin.peek()=='\n')cin.ignore();
		std::getline(cin,str);
		std::remove(str.begin(),str.end(),' ');
		cout<<"Heat? (y/n) ";
		heat=('y'==_get_restrict_che("yn"));
		cout<<endl;
		tube.assign(str,heat);
		tube.react_list(RXNS_FILE,obs_set);
		cout<<boost::lexical_cast<std::string>(tube)<<endl;
		tube.phys_obs(obs_set);
		for(it=obs_set.begin();it!=obs_set.end();++it)cout<<*it<<endl;
		obs_set.clear();
	}
	while(true);
}

void show_help()
{
	std::ifstream ifile(HELP_FILE);
	std::string str;
	std::getline(ifile,str,char(EOF));
	ifile.close();
	disp_msg(str);
}

void salt_analysis(char mode)
{
	std::string str;
	std::set<std::string> obs_set;
	std::set<std::string>::iterator it;
	tube_t tube;
	bool heat;
	comp_t comp;
	if(mode=='s')comp=get_random_salt();
	else if(mode=='a')comp=get_random_anion();
	else if(mode=='c')comp=get_random_cation();
	tube.add(comp);
	do
	{
		cout<<"\nEnter reagents:\t";
		if(cin.peek()=='\n')cin.ignore();
		getline(cin,str);
		if(str==comp.form)cout<<"Correct!"<<endl;
		else if(str=="answer")cout<<comp.form<<endl;
		else if(str=="new")
		{
			if(mode=='s')comp=get_random_salt();
			else if(mode=='a')comp=get_random_anion();
			else if(mode=='c')comp=get_random_cation();
			tube.comp_list.clear();
			tube.add(comp);
		}
		else if(str=="clear"){tube.comp_list.clear();tube.add(comp);}
		else if(str=="disp")cout<<boost::lexical_cast<std::string>(tube)<<endl;
		else if(str=="help")show_help();
		else
		{
			cout<<"Heat? (y/n) ";
			heat=('y'==_get_restrict_che("yn"));
			cout<<endl;
			tube.add(str,heat);
			tube.react_list(RXNS_FILE,obs_set);
			tube.phys_obs(obs_set);
			for(it=obs_set.begin();it!=obs_set.end();++it)cout<<*it<<endl;
			obs_set.clear();
		}
	}
	while(true);
}

int main()
{
	char ch;
	cout<<"Enter an option:\nr - Reactions\ns - Salt analysis\na - Anion analysis\nc - Cation analysis\n"<<endl;
	ch=_get_restrict_che("rsac");
	cout<<endl;
	if(ch=='r')reaction();
	else salt_analysis(ch);
	return 0;
}