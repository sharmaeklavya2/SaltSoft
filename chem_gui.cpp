#include<windows.h>
#include<string>
#include<commctrl.h>
#pragma comment(lib,"Comctl32")
#pragma comment(linker,"\"/manifestdependency:type='win32' \
name='Microsoft.Windows.Common-Controls' version='6.0.0.0' \
processorArchitecture='*' publicKeyToken='6595b64144ccf1df' language='*'\"")

void disp_error(const char* str)
{MessageBox(NULL,str,"Chem Analyze error",MB_ICONERROR|MB_OK);}

void disp_error(const std::string& str)
{MessageBox(NULL,str.c_str(),"Chem Analyze error",MB_ICONERROR|MB_OK);}

#include"chem.h"
#include"resource.h"

//#define DYNAMIC_REACT

void get_string_from_string_set(std::string& str,const std::set<std::string>& str_set)
{
	std::set<std::string>::const_iterator it=str_set.begin();
	str.clear();
	if(!str_set.empty()){str=*it;++it;}
	for(;it!=str_set.end();++it)str.append("\r\n").append(*it);
}

bool find_from_alias(const std::string& alias,std::string& actual)
{
	std::string str,str2;
	std::ifstream ifile(ALIAS_FILE);
	if(!ifile.is_open() || !ifile)
	{disp_error(std::string("Can't open ")+ALIAS_FILE);return 0;}
	do
	{
		std::getline(ifile,str,'\t');
		std::getline(ifile,str2,'\n');
		if(!ifile){ifile.close();return false;}
		if(str==alias){actual=str2;ifile.close();return true;}
	}
	while(true);
}

void get_string_from_edit(HWND hEdit,std::string& str)
{
	int size;
	char* buf;
	size=GetWindowTextLength(hEdit)+1;
	buf=new char[size];
	GetWindowText(hEdit,buf,size);
	str=buf;
	delete[] buf;
}

class listbox
{
	HWND hList;
public:
	listbox(HWND hlist=NULL):hList(hlist){}
	void set_handle(HWND hlist=NULL){hList=hlist;}
	int get_curr_index()const
	{return SendMessage(hList,LB_GETCURSEL,0,0);}
	void set_curr_index(int index)const
	{SendMessage(hList,LB_SETCURSEL,index,0);}
	void get_string(std::string& str,int index)
	{
		int size;
		char* buf;
		size=SendMessage(hList,LB_GETTEXTLEN,index,0)+1;
		if(size<=0)return;
		buf=new char[size];
		SendMessage(hList,LB_GETTEXT,index,(LPARAM)buf);
		str=buf;
		delete[] buf;
	}
	void delete_string(int index)const
	{SendMessage(hList,LB_DELETESTRING,index,0);}
	int add_string(const char* str)const
	{return SendMessage(hList,LB_ADDSTRING,0,(LPARAM)str);}
	int add_string(const std::string& str)const
	{return SendMessage(hList,LB_ADDSTRING,0,(LPARAM)str.c_str());}
	int add_data_from_file(const char* fname)const
	{
		std::string str;
		std::ifstream ifile(fname);
		int i=0;
		if(!ifile.is_open() || !ifile)
		{disp_error(std::string("Can't open ")+fname);return 0;}
		do
		{
			std::getline(ifile,str,'\t');
			ifile.ignore(1000,'\n');
			if(ifile){add_string(str);i++;}
			else break;
		}
		while(true);
		ifile.close();
		return i;
	}
	void add_comp_list(const std::set<comp_t>& comp_list)const
	{
		std::set<comp_t>::const_iterator it=comp_list.begin();
		for(;it!=comp_list.end();++it)
		{add_string(std::string(it->form)+"("+boost::lexical_cast<std::string>(it->st)+")");}
	}
	void clear()const
	{
		int noi=SendMessage(hList,LB_GETCOUNT,0,0);
		for(;noi>0;--noi)delete_string(0);
	}
};

BOOL CALLBACK DlgProc(HWND hDlg,UINT Msg,WPARAM wParam,LPARAM lParam)
{
	static listbox lbreg,lbreac;
	static HWND hEdit,hObs,hReact,hHeat;
	static tube_t tube;
	static std::string str,str2;
	int index;
	static std::set<std::string> obs_set;
	obs_set.clear();
	switch(Msg)
	{
	case WM_INITDIALOG:
		lbreg.set_handle(GetDlgItem(hDlg,IDC_REG));
		lbreac.set_handle(GetDlgItem(hDlg,IDC_REAC));
		hEdit=GetDlgItem(hDlg,IDC_EDIT);
		hObs=GetDlgItem(hDlg,IDC_POBS);
		hReact=GetDlgItem(hDlg,IDC_REACT);
		hHeat=GetDlgItem(hDlg,IDC_HEAT);
		lbreg.add_data_from_file(ALIAS_FILE);
		break;
	case WM_COMMAND:
		switch(LOWORD(wParam))
		{
		case IDCANCEL:
			EndDialog(hDlg,0);
			break;
		case IDC_ADD:
			get_string_from_edit(hEdit,str);
			tube.add(str,IsDlgButtonChecked(hDlg,IDC_HEAT)==BST_CHECKED);
#ifdef DYNAMIC_REACT
			SendMessage(hDlg,WM_COMMAND,IDC_REACT,0);
#else
			lbreac.clear();
			lbreac.add_comp_list(tube.comp_list);
#endif
			break;
		case IDC_CLEAR:
			tube.comp_list.clear();
			lbreac.clear();
			break;
		case IDC_REACT:
			tube.hight=(IsDlgButtonChecked(hDlg,IDC_HEAT)==BST_CHECKED);
			tube.react_list(RXNS_FILE,obs_set);
			tube.phys_obs(obs_set);
			get_string_from_string_set(str,obs_set);
			SendMessage(hObs,WM_SETTEXT,0,(LPARAM)str.c_str());
			lbreac.clear();
			lbreac.add_comp_list(tube.comp_list);
			break;
		case IDC_REMG:
			tube.remove_with_state(myenum::gas);
			lbreac.clear();
			lbreac.add_comp_list(tube.comp_list);
			break;
		case IDC_REML:
			tube.remove_with_state(myenum::liq);
			tube.remove_with_state(myenum::aq);
			lbreac.clear();
			lbreac.add_comp_list(tube.comp_list);
			break;
		case IDC_REMS:
			tube.remove_with_state(myenum::solid);
			lbreac.clear();
			lbreac.add_comp_list(tube.comp_list);
			break;
		case IDC_REG:
			switch(HIWORD(wParam))
			{
			case LBN_SELCHANGE:
				index=lbreg.get_curr_index();
				lbreg.get_string(str,index);
				find_from_alias(str,str2);
				SendMessage(hEdit,WM_SETTEXT,0,(LPARAM)str2.c_str());
				break;
			}
			break;
		}
		break;
	default:
		return FALSE;
	}
	return TRUE;
}

int __stdcall WinMain(HINSTANCE hInstance,HINSTANCE hPrevInstance,LPSTR lpCmdLine,int nShowCmd)
{
	DialogBox(hInstance,MAKEINTRESOURCE(IDD_RXN),NULL,DlgProc);
}