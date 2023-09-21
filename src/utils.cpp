#include"include/utils.h"
#include<vector>

using namespace std;
vector<std::vector<ListHead*>> grid;

void ListHead::addNew(Particle* par)
{
    //��������ͷ����һ�����Ӵ����Լ�����
    if (this->next == nullptr && this->tail != nullptr)
        cout << "something wrong\n";
    if (this->next == nullptr)
    {
        
        this->next = par;
        this->tail = par;
    }
    else
    {
        this->tail->next = par;
        par->last = this->tail;
        this->tail = this->tail->next;
    }
    this->num++;
}

ListHead::ListHead()
{
    this->next = nullptr;
    this->tail = nullptr;
    this->num = 0;
}
void initGrid(int num)
{
    for (int i = 0; i < num; i++)
    {
        vector<ListHead*> eachRow;
        for (int j = 0; j < num; j++)
        {
            auto block = new ListHead();
            eachRow.push_back(block);
        }
        grid.push_back(eachRow);
    }
}

void joinInGrid(Particle* par)
{
    //��һ�����Ӽ��뵽������
    if (par->i == -1 && par->j == -1)
    {
        //�޹��������ݵ�ǰ�����꽫����뵽ĳһ��Ͱ��
        par->i = (int)(par->x / 0.025) == 40 ? 39 : (int)(par->x / 0.025);
        par->j = (int)(par->y / 0.025) == 40 ? 39 : (int)(par->y / 0.025);
        //cout <<"no is: "<<debug++ << " par i is: " << par.i << " par j is: " << par.j << endl;
        grid[par->i][par->j]->addNew(par);
    }
    else
    {
        //�Ѿ�����ĳһ������
        int curr_i= (int)(par->x / 0.025) == 40 ? 39 : (int)(par->x / 0.025);
        int curr_j = (int)(par->y / 0.025) == 40 ? 39 : (int)(par->y / 0.025);
        if (curr_i == par->i && curr_j == par->j)
        {
            return;
        }
        else
        {
            //�����ӴӾɵ��������Ƴ��������뵽�µ�������
            removeFromOld(par);
            par->i = curr_i;
            par->j = curr_j;
            grid[par->i][par->j]->addNew(par);
        }
    }
    
}

void printGrid()
{
    for (int j = grid.size()-1; j >=0; j--)
    {
        for (int i = 0; i < grid.size(); i++)
        {
            cout << grid[i][j]->num << " ";
        }
        cout << endl;
    }
    cout << "**********************************************************************************************************************\n";
}

void removeFromOld(Particle* par)
{
    /*cout << "in function\n";
    cout << "i is: " << par->i << " and j is: " << par->j << endl;
    cout << "num is: " << grid[par->i][par->j]->num << endl;*/
    if (par->last == nullptr)
    {
        //�������ǵ�ǰ�����еĵ�һ������
        if (grid[par->i][par->j]->tail == par)
        {
            //��β��������
            grid[par->i][par->j]->tail = nullptr;
            grid[par->i][par->j]->next = nullptr;
            grid[par->i][par->j]->num--;
        }
        else
        {
            //���ײ����ǲ���β��������
            grid[par->i][par->j]->next = par->next;
            par->next->last = nullptr;
            par->next = nullptr;
            grid[par->i][par->j]->num--;
        }
    }
    else
    {
        //�����Ӳ��������еĵ�һ������
        if (grid[par->i][par->j]->tail == par)
        {
            //��β��������
            grid[par->i][par->j]->tail = par->last;
            par->last->next = nullptr;
            par->last = nullptr;
            par->next = nullptr;
            grid[par->i][par->j]->num--;
        }
        else
        {
            par->last->next = par->next;
            par->next->last = par->last;
            par->last = nullptr;
            par->next = nullptr;
            grid[par->i][par->j]->num--;
        }
    }
}

void checkFunc()
{
    initGrid(40);
    Particle* a = new Particle();
    a->x = 0.5;
    a->i = 0;
    a->j = 0;
    grid[0][0]->addNew(a);
    Particle* b = new Particle();
    b->x = 0.6;
    b->i = 0;
    b->j = 0;
    grid[0][0]->addNew(b);
    printGrid();
    auto itea = grid[0][0]->next;
    while (itea != nullptr)
    {
        cout << itea->x << endl;
        itea = itea->next;
    }
    removeFromOld(b);
    printGrid();
}