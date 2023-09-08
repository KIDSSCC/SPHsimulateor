#include"include/utils.h"
#include<vector>

using namespace std;
vector<std::vector<ListHead*>> grid;

void ListHead::addNew(Particle* par)
{
    //对于链表头，将一个粒子串到自己身上
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
    //将一个粒子加入到网格中
    if (par->i == -1 && par->j == -1)
    {
        //无归属，根据当前的坐标将其加入到某一个桶中
        par->i = (int)(par->x / 0.025) == 40 ? 39 : (int)(par->x / 0.025);
        par->j = (int)(par->y / 0.025) == 40 ? 39 : (int)(par->y / 0.025);
        //cout <<"no is: "<<debug++ << " par i is: " << par.i << " par j is: " << par.j << endl;
        grid[par->i][par->j]->addNew(par);
    }
    else
    {
        //已经属于某一个网格
        int curr_i= (int)(par->x / 0.025) == 40 ? 39 : (int)(par->x / 0.025);
        int curr_j = (int)(par->y / 0.025) == 40 ? 39 : (int)(par->y / 0.025);
        if (curr_i == par->i && curr_j == par->j)
        {
            return;
        }
        else
        {
            //将粒子从旧的网格中移出来，加入到新的网格中
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
        //该粒子是当前网格中的第一个粒子
        if (grid[par->i][par->j]->tail == par)
        {
            //是尾部的粒子
            grid[par->i][par->j]->tail = nullptr;
            grid[par->i][par->j]->next = nullptr;
            grid[par->i][par->j]->num--;
        }
        else
        {
            //是首部但是不是尾部的粒子
            grid[par->i][par->j]->next = par->next;
            par->next->last = nullptr;
            par->next = nullptr;
            grid[par->i][par->j]->num--;
        }
    }
    else
    {
        //该粒子不是网格中的第一个粒子
        if (grid[par->i][par->j]->tail == par)
        {
            //是尾部的粒子
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