/*This is a program to simulate chemical tests. A list of reactions will be provided (inbuilt[+user]).
The user will enter products and reactants will be computed.
Physical characteristics of the compound will also be computed using files
Files used by the program are in a folder called progdata

When I started making this project, I didn't know that usually STL uses
(!(a>b) && !(b<a)) instead of a==b
Because comp_t::operator==(const comp_t&)const is not a transitive function,
I had to write my own algorithms at some places.
Also, I cannot always use the fact that std::set<comp_t> is sorted
as MgCO3(s)==CO32-(aq) while sorting is on the basis of chemical formula.
Hence, while checking if one set includes the other,
I can't use the algorithm optimized for sorted ranges

Because of the above reason, the program may seem a bit chaotic.
*/

#pragma warning(disable:4996)

#include<iostream>
#include<fstream>
#include<cstring>
#include<cstdlib>
#include<ctime>
#include<conio.h>	//for _get_restrict_che()
#include<string>
#include<set>
#include<algorithm>
#include<boost\foreach.hpp>
#include<boost\tokenizer.hpp>
#include<boost\lexical_cast.hpp>
//#include<windows.h>

#define COMP_FORM_SIZE 20
#define FILE_PATH_SIZE 256
#define FILE_NAME_SIZE 30

//Std File Names
#define CATION_LIST "progdata\\cations.txt"
#define ANION_LIST	"progdata\\anions.txt"
#define SALT_SPLIT	"progdata\\salt split.txt"
#define PPTS	"progdata\\precipitates.txt"
//#define PPTS_FOLDER	"progdata\\ppt\\";
#define NEVER_PPT_ION	"progdata\\never ppt ion.txt"
#define DISSOC_MEDIA	"progdata\\dissoc media.txt"
#define RXNS_FOLDER		"progdata\\rxns\\"
#define RXNS_FILE		"progdata\\rxns.txt"
#define ION_COLOR	"progdata\\coloured ions.txt"
#define GAS_SMELL_COLOR	"progdata\\gases.txt"
#define ANAL_SALTS	"progdata\\salt analysis\\salts.txt"
#define ANAL_CATIONS	"progdata\\salt analysis\\cations.txt"
#define ANAL_ANIONS	"progdata\\salt analysis\\anions.txt"
#define HELP_FILE	"progdata\\help.txt"

//global=functions================================================================

int crude_split_salt(const char* salt,char* cation,char* anion)
/*splits salt into cation and anion
the cation can't be polyatomic
the anion isn't correctly identified.
Eg. AlCl3 as salt gives Cl3 as anion instead of Cl
Because of the above constraints, the function is named 'crude'
The cations and anions are actually found by comp_t::split()
comp_t::split() calls this function as a subprocedure*/
{
	int i,j,k,n=std::strlen(salt),ret;
	char num[10];
	bool cbrac_present=false;
	cation[0]='\0';anion[0]='\0';
	
	if(salt[0]<'A' || salt[0]>'Z')return 0;
	for(i=1;i<n && salt[i]>='a' && salt[i]<='z';i++);
	//keep looping till an non-lowercase character is encountered
	for(j=0;j<i;j++)cation[j]=salt[j];cation[j]='\0';
	//keep the first i characters in cation
	for(j=0;i<n && salt[i]>='0' && salt[i]<='9';i++,j++)num[j]=salt[i];
	num[j]='\0';
	if(strlen(num)==0)ret=1;
	else ret=std::atoi(num);
	//count number of cations

	if(salt[i]=='(')
	{j=i+1;for(i=n-1;i>=0 && salt[i]!=')';--i);}
	else
	{j=i;/*for(;i>=0 && salt[i]>='0' && salt[i]<='9';--i);i++;*/i=n;}
	for(k=0;j<i;j++,k++)anion[k]=salt[j];
	anion[k]='\0';
	return ret;
}

bool ion_match(const char* ion,const char* stdion)
//checks if anion obtained after crude_salt_split is correct
{
	int i,h=strlen(ion),hstd=strlen(stdion);
	if(h<hstd)return false;
	for(i=0;i<hstd;i++)if(ion[i]!=stdion[i])return false;
	return true;
}

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

std::istream& my_getline(std::istream& is,std::string& str,char delim='\n')
{
	if(is.peek()=='\n')is.ignore();
	std::getline(is,str,delim);
	return is;
}
std::istream& my_getline(std::istream& is,char* str,int size,char delim='\n')
{
	if(is.peek()=='\n')is.ignore();
	is.getline(str,size,delim);
	return is;
}
//Classes=========================================================================

namespace myenum
{
	enum state_t{gas,solid,liq,aq,dil,conc,vdil,nil,any_st};
	//enum containing possible states of compound
	//enum charge_t{positive,negative,neutral};
	//enum containing possible charges on a species
	enum comp_type_t{salt,cation,anion,notsalt,none};
	//Compound type:	electrovalent, covalent, unspecified
}

class comp_t;
class tube_t;
class rxn_t;

namespace boost
{
	template<>
	myenum::state_t lexical_cast<myenum::state_t,std::string>(const std::string& str)
	{
		if(str=="g")return myenum::gas;
		else if(str=="l")return myenum::liq;
		else if(str=="s")return myenum::solid;
		else if(str=="aq")return myenum::aq;
		else if(str=="dil")return myenum::dil;
		else if(str=="conc")return myenum::conc;
		else if(str=="vdil")return myenum::vdil;
		else if(str=="any")return myenum::any_st;
		else return myenum::nil;
	}
	template<>
	std::string lexical_cast<std::string,myenum::state_t>(const myenum::state_t& state)
	{
		if(state==myenum::solid)return "s";
		else if(state==myenum::aq)return "aq";
		else if(state==myenum::liq)return "l";
		else if(state==myenum::gas)return "g";
		else if(state==myenum::dil)return "dil";
		else if(state==myenum::conc)return "conc";
		else if(state==myenum::vdil)return "vdil";
		else if(state==myenum::any_st)return "any";
		else return "nil";
	}
	template<>
	myenum::comp_type_t lexical_cast<myenum::comp_type_t,std::string>(const std::string& str)
	{
		if(str=="salt")return myenum::salt;
		else if(str=="notsalt" || str=="not salt")return myenum::notsalt;
		else if(str=="cation")return myenum::cation;
		else if(str=="anion")return myenum::anion;
		else return myenum::none;
	}
	template<>
	std::string lexical_cast<std::string,myenum::comp_type_t>(const myenum::comp_type_t& comp_type)
	{
		if(comp_type==myenum::salt)return "salt";
		else if(comp_type==myenum::notsalt)return "notsalt";
		else if(comp_type==myenum::cation)return "cation";
		else if(comp_type==myenum::anion)return "anion";
		else return "none";
	}
	template<>
	comp_t lexical_cast<comp_t,std::string>(const std::string& form_with_state);
	template<>
	std::string lexical_cast<std::string,comp_t>(const comp_t& comp);
	template<>
	std::string lexical_cast<std::string,tube_t>(const tube_t& comp);
	template<>
	rxn_t lexical_cast<rxn_t,std::string>(const std::string& str);
}

//--------------------------------------------------------------------------------

class comp_t
//compound class
//contains formula and state of compound
{public:
	char form[COMP_FORM_SIZE];	//chemical formula of compound object
	char cation[COMP_FORM_SIZE];	//cation of compound (if exists)
	char anion[COMP_FORM_SIZE];		//anion of compound (if exists)
	myenum::comp_type_t comp_type;
	myenum::state_t st;	//physical state of compound
	comp_t():st(myenum::nil),comp_type(myenum::none)
	{form[0]='\0';cation[0]='\0';anion[0]='\0';}
	comp_t(const char* formula,myenum::state_t state,myenum::comp_type_t type=myenum::none):st(state),comp_type(type)
	{std::strcpy(form,formula);cation[0]='\0';anion[0]='\0';}
	comp_t(const comp_t& comp):st(comp.st),comp_type(comp.comp_type)
	{
		std::strcpy(form,comp.form);
		std::strcpy(cation,comp.cation);
		std::strcpy(anion,comp.anion);
	}
	bool operator<(const comp_t& comp)const;
	//used for compatibility with std::set<comp_t>
	bool operator==(const comp_t& comp)const
	{
		return (std::strcmp(form,comp.form)==0
			|| (anion[0]!='\0' && std::strcmp(anion,comp.form)==0)
			|| (cation[0]!='\0' && std::strcmp(cation,comp.form)==0))
			&& (st==comp.st || st==myenum::any_st || comp.st==myenum::any_st);
	}
	bool exact_equal(const comp_t& comp)const
	{
		return std::strcmp(form,comp.form)==0
			&& (st==comp.st || st==myenum::any_st || comp.st==myenum::any_st);
	}
	bool ion_equal(const comp_t& comp)const
	{
		return ((anion[0]!='\0' && std::strcmp(anion,comp.form)==0)
			|| (cation[0]!='\0' && std::strcmp(cation,comp.form)==0))
			&& (st==comp.st || st==myenum::any_st || comp.st==myenum::any_st);
	}
	bool anion_equal(const comp_t& comp)const
	{
		return anion[0]!='\0' && std::strcmp(anion,comp.form)==0
			&& (st==comp.st || st==myenum::any_st || comp.st==myenum::any_st);
	}
	bool cation_equal(const comp_t& comp)const
	{
		return cation[0]!='\0' && std::strcmp(cation,comp.form)==0
			&& (st==comp.st || st==myenum::any_st || comp.st==myenum::any_st);
	}
	void operator=(const std::string& form_with_state);
	//eg. input "(NH4)2CO3(s)" will set form as "(NH4)2CO3" and state as solid.
	void split();
	//gets the cation and anion of salt
	//leaves them "" if compound is not a salt
	//sets comp_type as myenum::salt if salt and myenum::notsalt if not a salt
};

bool comp_t::operator<(const comp_t& comp)const
{
	int streq=strcmp(form,comp.form);
	if(streq<0)return true;
	else if(streq>0)return false;
	else return st<comp.st;
}

//--------------------------------------------------------------------------------

class tube_t
//test tube class
//contains list of compounds in it and temperature
{public:
	std::set<comp_t> comp_list;	//list of compounds in test tube
	bool hight;	//is temperature high?
	bool show_mech;	//
	tube_t(bool high_temp=false,bool show_mechanism=false):
	hight(high_temp),show_mech(show_mechanism){}
	void add(const comp_t& comp)
	{comp_list.insert(comp);}
	void add(const char* formula,myenum::state_t state,myenum::comp_type_t type=myenum::none)
	{comp_list.insert(comp_t(formula,state,type));}
	bool carry_rxn(const rxn_t& rxn);
	//carries out given reaction in tube (if possible)
	//rxn should be pre-filled using rxn_t::operator=(const std::string& str);
	//return true if reaction happened
	void react_prep();
	//dissociate double salts and other reaction initialization
	void assign(std::string str,bool heat=false);
	void add(std::string str,bool heat=false);
	void precipitate();
	//checks for ions present and precipitates precipitable salts.
	void react(const char* fname);
	//carries out reactions from a reaction file
	void react_list(const char* fname=RXNS_FILE);
	//carries out reactions from a reaction list file
	//The reaction list file contains file names from the RXNS_FOLDER
	//These file names have reactions in them
	//displays computer mechanism of reaction using std::cout if show_mech is true
	void phys_obs_disp()const;
	//physical observation of test tube
	//disp_each displays each observation separately using disp_msg
};

//--------------------------------------------------------------------------------

class rxn_t
//reaction class
//contains list of compounds in it and temperature
{public:
	std::set<comp_t> reac_list;	//list of reactants used
	std::set<comp_t> prod_list;	//list of products formed
	std::set<comp_t> cat_list;	//list of catalysts
	std::string comment;
	bool hight;	//is temperature high
	bool lowt;	//is temperature low
	rxn_t(bool high_temp=false,bool low_temp=false):hight(high_temp),lowt(low_temp){}
	void operator=(const std::string& str);
	//BUG: interprets '+' sign in anions as '+' of adding reagents
};

//Implementation==================================================================

namespace boost
{
	template<>
	comp_t lexical_cast<comp_t,std::string>(const std::string& form_with_state)
	{
		int pos;std::string state;
		comp_t comp;
		//find last occurence of '('
		comp.comp_type=myenum::none;
		pos=form_with_state.find_last_of("(");
		if(pos==-1)
		{comp.st=myenum::any_st;std::strcpy(comp.form,form_with_state.c_str());}
		else 
		{
			state.assign(form_with_state,pos+1,form_with_state.size()-pos-2);
			comp.st=boost::lexical_cast<myenum::state_t>(state);
			if(comp.st==myenum::nil)
			{comp.st=myenum::any_st;std::strcpy(comp.form,form_with_state.c_str());}
			else std::strcpy(comp.form,form_with_state.substr(0,pos).c_str());
		}
		//assign charge
		if(comp.form[pos-1]=='+')comp.comp_type=myenum::cation;
		else if(comp.form[pos-1]=='-')comp.comp_type=myenum::anion;
		else comp.comp_type=myenum::none;
		//remove coefficients
		for(pos=0;comp.form[pos]>='0' && comp.form[pos]<='9';++pos);
		std::strcpy(comp.form,comp.form+pos);
		return comp;
	}
	template<>
	std::string lexical_cast<std::string,comp_t>(const comp_t& comp)
	{
		std::string str(comp.form);
		str.append(1,'(');
		str.append(boost::lexical_cast<std::string>(comp.st));
		str.append(1,')');
		return str;
	}
	template<>
	std::string lexical_cast<std::string,tube_t>(const tube_t& tube)
	{
		std::string str;
		if(tube.comp_list.empty())return "";
		std::set<comp_t>::const_iterator it=tube.comp_list.begin();
		for(;it!=tube.comp_list.end();++it)
		{
			str.append(lexical_cast<std::string>(*it));
			str.append(1,'+');
		}
		str.erase(str.size()-1,1);
		return str;
	}
	template<>
	rxn_t lexical_cast<rxn_t,std::string>(const std::string& str2)
	{
		int pos1,pos2;
		std::string reac,prod,cat,temp,str=str2;
		comp_t comp;
		rxn_t rxn;
		int i;
		boost::char_separator<char> sep("&");
		boost::tokenizer<boost::char_separator<char> >tokens(reac);

		//take out comment
		pos1=str.find("//");
		if(pos1==-1)pos1=str.find("\\\\");
		else if(pos1!=-1)
		{
			rxn.comment.assign(str,pos1+2,std::string::npos);
			str.erase(pos1,std::string::npos);
			for(i=0;i<int(rxn.comment.size());i++)
			{if(rxn.comment[i]=='_')rxn.comment[i]=' ';}
		}
		//splits reaction into reactants,products,catalysts
		pos1=str.find("--");
		reac.assign(str,0,pos1);
		pos2=str.find("-->",pos1);
		if(pos2-pos1>2)cat.assign(str,pos1+2,pos2-pos1-2);
		prod.assign(str,pos2+3,std::string::npos);

		//split reactants
		if(reac.size()>0)
		{
			for(i=0;i<int(reac.size())-1;i++)if(reac[i]=='+' && ((i>0 && reac[i-1]==')') || (reac[i+1]!='+' && reac[i+1]!='(')))reac[i]='&';
			tokens.assign(reac,sep);
			BOOST_FOREACH(temp,tokens)
			{rxn.reac_list.insert(boost::lexical_cast<comp_t>(temp));}
		}

		//split products
		if(prod.size()>0)
		{
			for(i=0;i<int(prod.size())-1;i++)if(prod[i]=='+' && ((i>0 && prod[i-1]==')') || (prod[i+1]!='+' && prod[i+1]!='(')))prod[i]='&';
			tokens.assign(prod,sep);
			BOOST_FOREACH(temp,tokens)
			{rxn.prod_list.insert(boost::lexical_cast<comp_t>(temp));}
		}

		//split catalysts (if exist)
		if(pos2-pos1>2)
		{
			for(i=0;i<int(cat.size())-1;i++)if(cat[i]=='+' && ((i>0 && cat[i-1]==')') || (cat[i+1]!='+' && cat[i+1]!='(')))cat[i]='&';
			tokens.assign(cat,sep);
			BOOST_FOREACH(temp,tokens)
			{
				if(temp=="heat"){rxn.hight=true;rxn.lowt=false;}
				else if(temp=="cold"){rxn.lowt=true;rxn.hight=false;}
				else rxn.cat_list.insert(boost::lexical_cast<comp_t>(temp));
			}
		}
		return rxn;
	}
}

void get_random_line(const char* fname,char* line)
//gets a random line from file 'fname' into 'line'
{
	std::ifstream ifile(fname);
	int i,size,pos;
	comp_t comp;
	if(!ifile.is_open())disp_error(std::string("Can't open ")+fname);
	for(size=0;ifile.ignore(10000,'\n');++size);
	ifile.clear();
	ifile.seekg(0);
	srand(unsigned(std::time(NULL)));pos=rand()%size;
	for(i=0;i<=pos;i++)ifile>>line;
	ifile.close();
}

bool my_includes(const std::set<comp_t>& bigset,const std::set<comp_t>& smallset)
/*my own version of std::includes written specially for std::set<comp_t>
This had to be done because comp_t::operator==(const comp_t&)const was not transitive
So I couldn't rewrite comp_t::operator<(const comp_t&)const so that STL can
use (!(a<b) && !(b<a)) to get the same result as a==b*/
{
	std::set<comp_t>::const_iterator its,itb;
	bool found_elem_in_big=false;
	if(smallset.size()>bigset.size())return false;
	for(its=smallset.begin();its!=smallset.end();++its)
	{
		for(itb=bigset.begin();itb!=bigset.end();++itb)
		{
			found_elem_in_big=itb->operator==(*its);
			if(found_elem_in_big)break;
		}
		if(!found_elem_in_big)return false;
	}
	return true;
}

void comp_t::operator=(const std::string& form_with_state)
{
	int pos;std::string state;
	//find last occurence of '('
	comp_type=myenum::none;
	pos=form_with_state.find_last_of("(");
	if(pos==-1)
	{st=myenum::any_st;std::strcpy(form,form_with_state.c_str());}
	else 
	{
		state.assign(form_with_state,pos+1,form_with_state.size()-pos-2);
		st=boost::lexical_cast<myenum::state_t>(state);
		if(st==myenum::nil)
		{st=myenum::any_st;std::strcpy(form,form_with_state.c_str());}
		else std::strcpy(form,form_with_state.substr(0,pos).c_str());
	}
	//assign charge
	if(form[pos-1]=='+')comp_type=myenum::cation;
	else if(form[pos-1]=='-')comp_type=myenum::anion;
	else comp_type=myenum::none;
	//remove coefficients
	for(pos=0;form[pos]>='0' && form[pos]<='9';++pos);
	std::strcpy(form,form+pos);
}

void comp_t::split()
{
	char stdform[COMP_FORM_SIZE],stdcat[COMP_FORM_SIZE],stdan[COMP_FORM_SIZE],charge[4];
	bool found;
	//first check "salt split.txt" for standard splitting
	std::ifstream ifile(SALT_SPLIT);
	if(!ifile.is_open())disp_error(std::string("Can't open")+SALT_SPLIT);
	while(ifile>>stdform>>stdcat>>stdan)if(std::strcmp(stdform,form)==0)
	{
		ifile.close();
		std::strcpy(cation,stdcat);
		std::strcpy(anion,stdan);
		comp_type=myenum::salt;
		return;
	}
	ifile.close();
	//if not found, use call crude_split_salt(form,cation,anion) to get crude ions
	crude_split_salt(form,cation,anion);
	//then check if cation exists using "cations.txt"
	//if it exists, get its charge and append the charge at the end of the data member 'cation'
	ifile.clear();
	ifile.open(CATION_LIST);
	if(!ifile.is_open())disp_error(std::string("Can't open")+CATION_LIST);
	found=false;
	while(ifile>>stdcat>>charge)if(std::strcmp(stdcat,cation)==0)
	{
		ifile.close();
		found=true;
		std::strcat(cation,charge);
		break;
	}
	ifile.close();
	if(!found){cation[0]='\0';anion[0]='\0';comp_type=myenum::notsalt;return;}
	//check if anion exists using "anions.txt" and function ion_match()
	//if it exists, get its charge and append the charge at the end of the data member 'anion'
	ifile.clear();
	ifile.open(ANION_LIST);
	if(!ifile.is_open())disp_error(std::string("Can't open")+ANION_LIST);
	found=false;
	while(ifile>>stdan>>charge)if(ion_match(anion,stdan))
	{
		ifile.close();
		found=true;
		std::strcpy(anion,stdan);
		std::strcat(anion,charge);
		break;
	}
	ifile.close();
	if(!found)
	{
		cation[0]='\0';
		anion[0]='\0';
		comp_type=myenum::notsalt;
		return;
	}
	comp_type=myenum::salt;
}
bool tube_t::carry_rxn(const rxn_t& rxn)
{
	std::set<comp_t>::const_iterator itr;
	std::set<comp_t>::iterator it;
	comp_t comp;
	if(!my_includes(comp_list,rxn.reac_list))return false;
	if(!my_includes(comp_list,rxn.cat_list))return false;
	if((rxn.hight && !hight)||(rxn.lowt && hight))return false;
	for(itr=rxn.reac_list.begin();itr!=rxn.reac_list.end();++itr)
	{
		for(it=comp_list.begin();it!=comp_list.end();++it)
		{
			if(it->exact_equal(*itr))
			{comp_list.erase(it);break;}
			else if(it->anion_equal(*itr))
			{
				comp=comp_t(it->cation,myenum::aq,myenum::cation);
				comp_list.erase(it);
				comp_list.insert(comp);
				break;
			}
			else if(it->cation_equal(*itr))
			{
				comp=comp_t(it->anion,myenum::aq,myenum::anion);
				comp_list.erase(it);
				comp_list.insert(comp);
				break;
			}
		}
	}
	comp_list.insert(rxn.prod_list.begin(),rxn.prod_list.end());
	if(!rxn.comment.empty())
	{disp_msg(rxn.comment);}
	if(show_mech)std::cout<<boost::lexical_cast<std::string>(*this)<<std::endl;
	return true;
}

void tube_t::assign(std::string str,bool heat/*=false*/)
{
	std::string comp_str;
	comp_list.clear();
	hight=heat;
	int i;
	for(i=0;i<int(str.size())-1;i++)if(str[i]=='+' && ((i>0 && str[i-1]==')') || (str[i+1]!='+' && str[i+1]!='(')))str[i]='&';
	boost::char_separator<char> sep("&");
	boost::tokenizer<boost::char_separator<char> >tokens(str,sep);
	BOOST_FOREACH(comp_str,tokens)
	{comp_list.insert(boost::lexical_cast<comp_t>(comp_str));}
}

void tube_t::add(std::string str,bool heat/*=false*/)
{
	std::string comp_str;
	hight=heat;
	int i;
	for(i=0;i<int(str.size())-1;i++)if(str[i]=='+' && ((i>0 && str[i-1]==')') || (str[i+1]!='+' && str[i+1]!='(')))str[i]='&';
	boost::char_separator<char> sep("&");
	boost::tokenizer<boost::char_separator<char> >tokens(str,sep);
	BOOST_FOREACH(comp_str,tokens)
	{comp_list.insert(boost::lexical_cast<comp_t>(comp_str));}
}

void tube_t::react_prep()
{
	char salt[COMP_FORM_SIZE],dummy[COMP_FORM_SIZE];
	bool is_ppt=false,found=false,never_ppt=false;
	comp_t comp;
	std::set<comp_t> new_comp_list;
	std::ifstream ifile,ifile2;
	std::set<comp_t>::iterator it;
	char ch;
	//find type of and cations and anions in each compound
	for(it=comp_list.begin();it!=comp_list.end();++it)
	{
		if(it->comp_type==myenum::none)
		{it->split();}
	}
	//check for dissoc media
	ifile.open(DISSOC_MEDIA);
	if(!ifile.is_open())disp_error(std::string("Can't open")+DISSOC_MEDIA);
	while(ifile>>dummy)
	{
		comp=dummy;
		found=std::binary_search(comp_list.begin(),comp_list.end(),comp);
		if(found)break;
	}
	ifile.close();
	if(!found)return;
	//convert solid ions to aqueous ions
	for(it=comp_list.begin();it!=comp_list.end();++it)
	{if(it->comp_type==myenum::anion || it->comp_type==myenum::cation)it->st=myenum::aq;}
	//dissociate if not in ppts file and found in NEVER_PPT_ION
	ifile.clear();
	ifile.open(PPTS);
	ifile2.open(NEVER_PPT_ION);
	if(!ifile.is_open())disp_error(std::string("Can't open")+PPTS);
	if(!ifile2.is_open())disp_error(std::string("Can't open")+NEVER_PPT_ION);
	new_comp_list=comp_list;
	BOOST_FOREACH(comp,comp_list)if(comp.comp_type==myenum::salt)
	{
		is_ppt=false;never_ppt=false;
		while(ifile2>>dummy)
		{
			ch=dummy[strlen(dummy)-1];
			if((ch=='+' && strcmp(comp.cation,dummy)==0) || (ch=='-' && strcmp(comp.anion,dummy)==0)) 
			{never_ppt=true;break;}
		}
		if(!never_ppt)while(ifile>>salt>>dummy)
		{
			is_ppt=(std::strcmp(salt,comp.form)==0);
			if(is_ppt)break;
		}
		if(!is_ppt)
		{
			new_comp_list.insert(comp_t(comp.cation,myenum::aq,myenum::cation));
			new_comp_list.insert(comp_t(comp.anion,myenum::aq,myenum::anion));
			new_comp_list.erase(comp);
		}
	}
	comp_list=new_comp_list;
	ifile.close();
}

void tube_t::precipitate()
{
	comp_t comp,compcat,compan;
	std::string ppt,dummy;
	std::ifstream ifile(PPTS);
	std::set<comp_t>::iterator itcat,itan;
	bool fcat,fan;
	while(ifile>>ppt>>dummy)
	{
		comp=ppt;
		comp.st=myenum::solid;
		comp.split();
		compcat=comp_t(comp.cation,myenum::aq,myenum::cation);
		compan=comp_t(comp.anion,myenum::aq,myenum::anion);
		itcat=comp_list.find(compcat);
		fcat=itcat!=comp_list.end();
		itan=comp_list.find(compan);
		fan=itan!=comp_list.end();
		if(fcat&&fan)
		{
			comp_list.erase(compcat);
			comp_list.erase(compan);
			comp_list.insert(comp);
		}
	}
}


/*void tube_t::ppt_by_cation()
{
	std::ifstream ifile;
	char path[FILE_PATH_SIZE],ion[FILE_NAME_SIZE];
	std::set<comp_t>::iterator it,it2;
	std::set<comp_t>::to_remove,to_add;
	for(it=comp_list.begin();it!=comp_list.end();++it)
	{
		if(!it->comp_type==myenum::cation)continue;
		std::strcpy(ion,it->form);
		ifile.clear();
		std::strcpy(path,PPTS_FOLDER);
		std::strcat(path,ion)
		ifile.open(path);
		if(!ifile.open()){ifile.close();continue;}
		to_remove.insert(*it);
		while(ifile>>ion)
		{
			for(it2=comp_list.begin();it2!=comp_list.end();++it2)
			{if(it2->comp_type==myenum::anion && std::strcmp(it2->form,ion)==0)to_remove.insert(*it2);}
			//generate ppt from *it and *it2 and insert it in 'to_add'
		}
		ifile.close();
	}
	//remove 'to_remove' from 'comp_list'
	//add 'to_add' to 'comp_list'
}*/

void tube_t::react(const char* fname)
{
	std::ifstream ifile(fname);
	if(!ifile.is_open())disp_error(std::string("Can't open ")+fname);
	std::string str;
	while(ifile>>str)
	{carry_rxn(boost::lexical_cast<rxn_t>(str));}
	ifile.close();
}

void tube_t::react_list(const char* fname/*=RXNS_FILE*/)
{
	char path[FILE_PATH_SIZE],str[FILE_NAME_SIZE];
	std::ifstream ifile(fname);
	if(!ifile.is_open())disp_error(std::string("Can't open ")+fname);
	react_prep();
	while(my_getline(ifile,str,FILE_NAME_SIZE))
	{
		strcpy(path,RXNS_FOLDER);
		strcat(path,str);
		react(path);
	}
	precipitate();
}

void tube_t::phys_obs_disp()const
{
	std::ifstream pptfile(PPTS);
	std::ifstream ionfile(ION_COLOR);
	std::ifstream gasfile(GAS_SMELL_COLOR);
	std::set<comp_t>::const_iterator it=comp_list.begin();
	std::string comp,color,smell;
	for(;it!=comp_list.end();++it)
	{
		if(it->st==myenum::solid)
		{
			pptfile.clear();
			pptfile.seekg(0);
			while(pptfile>>comp>>color)if(comp==it->form)disp_msg(color+" ppt");
		}
		else if(it->st==myenum::gas)
		{
			gasfile.clear();
			gasfile.seekg(0);
			while(gasfile>>comp>>smell>>color)if(comp==it->form)
			{
				if(color!="none")disp_msg(color+" gas");
				if(smell!="none")disp_msg(smell+" smell");
			}
		}
		else if(it->comp_type==myenum::anion || it->comp_type==myenum::cation || it->comp_type==myenum::notsalt)
		{
			ionfile.clear();
			ionfile.seekg(0);
			while(ionfile>>comp>>color)if(comp==it->form)disp_msg(color+" solution");
		}
	}
	pptfile.close();
	ionfile.close();
	gasfile.close();
}

//================================================================================

using namespace std;

comp_t get_random_anion()
{
	comp_t comp;
	get_random_line(ANAL_ANIONS,comp.form);
	comp.st=myenum::solid;
	comp.comp_type=myenum::anion;
	return comp;
}

comp_t get_random_cation()
{
	comp_t comp;
	get_random_line(ANAL_CATIONS,comp.form);
	comp.st=myenum::solid;
	comp.comp_type=myenum::cation;
	return comp;
}

comp_t get_random_salt()
{
	comp_t comp;
	get_random_line(ANAL_SALTS,comp.form);
	comp.st=myenum::solid;
	comp.comp_type=myenum::salt;
	comp.split();
	return comp;
}

void reaction()
{
	std::string str;
	tube_t tube;
	bool heat;
	do
	{
		cout<<"\nEnter reactants:\t";
		cin>>str;
		cout<<"Heat? (y/n) ";
		heat=('y'==_get_restrict_che("yn"));
		cout<<endl;
		tube.assign(str,heat);
		tube.react_list(RXNS_FILE);
		cout<<boost::lexical_cast<std::string>(tube)<<endl;
		tube.phys_obs_disp();
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
			tube.react_list(RXNS_FILE);
			tube.phys_obs_disp();
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