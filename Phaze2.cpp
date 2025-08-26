#include <iostream>
#include <vector>
#include <SDL2/SDL.h>
#include <SDL2/SDL2_gfx.h>
#include <SDL2/SDL_ttf.h>
#include <SDL2/SDL_image.h>
#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <cstdlib>
#include <clocale>
#include <windows.h>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <cstdlib>   // strtod
#include <clocale>   // setlocale
#include "phase1_api.h"
#include <shellapi.h> // اضافه شود در بالای فایل
#include <commdlg.h>  // برای GetSaveFileNameA و فلگ‌های مربوطه
#include <commdlg.h>  // دیالوگ ذخیرهٔ ویندوز (GetSaveFileNameA)
#include <sstream>
#include <cstdio>      // freopen_s, fflush
#include <io.h>        // _dup, _dup2, _fileno   (ویندوز)
#include <fcntl.h>     // _O_* flags (اگر لازم شد)
#include <regex>

// برمی‌گرداند اندیسی که data[idx].first <= x < data[idx+1].first
// اگر x بیرون بازه بود، در لبه‌ها clamp می‌کند.
static size_t bracket_index(const std::vector<std::pair<double,double>>& data, double x) {
    if (data.size() <= 1) return 0;
    if (x <= data.front().first) return 0;
    if (x >= data.back().first)  return data.size() - 2;
    size_t lo = 0, hi = data.size() - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) >> 1;
        if (data[mid].first <= x) lo = mid; else hi = mid;
    }
    return lo;
}
static std::string elName(size_t idx) ;

static int pickNodeAt(int mx, int my) ;
struct CursorUI {
    bool enabled = false;   // حالت نشانگر فعال/غیرفعال
    bool hasPoint = false;  // آیا نقطه‌ای انتخاب شده؟
    int  px = 0, py = 0;    // مختصات پیکسلی نقطه انتخاب شده
    size_t idx = 0;         // اندیس داده انتخاب شده
    SDL_Rect icon{8, 8, 18, 18}; // آیکون کوچک گوشه بالا-چپ
};

static void draw_cursor_icon(SDL_Renderer* r, const SDL_Rect& rc, bool on) {
    // قاب
    SDL_SetRenderDrawColor(r, on ? 0 : 60, on ? 160 : 60, on ? 0 : 60, 255);
    SDL_RenderDrawRect(r, &rc);
    // علامت +
    int cx = rc.x + rc.w/2, cy = rc.y + rc.h/2;
    SDL_RenderDrawLine(r, cx-5, cy,   cx+5, cy);
    SDL_RenderDrawLine(r, cx,   cy-5, cx,   cy+5);
}
static void cleanupDanglingNodes() ;
static bool point_in_rect(int x, int y, const SDL_Rect& rc) {
    return x >= rc.x && x <= rc.x + rc.w && y >= rc.y && y <= rc.y + rc.h;
}

static void draw_grid_linearXY(SDL_Renderer* r,
                               int L, int T, int W, int H, int R, int B,
                               int nx = 20, int ny = 20) {
    SDL_SetRenderDrawColor(r, 220,220,220,255); // خاکستری روشن
    // محدودهٔ داخل نمودار
    int x0 = L, x1 = W - R;
    int y0 = T, y1 = H - B;

    // خطوط عمودی
    for (int i = 1; i < nx; ++i) {
        int x = x0 + (x1 - x0) * i / nx;
        SDL_RenderDrawLine(r, x, y0, x, y1);
    }
    // خطوط افقی
    for (int j = 1; j < ny; ++j) {
        int y = y1 - (y1 - y0) * j / ny;
        SDL_RenderDrawLine(r, x0, y, x1, y);
    }
}
static void draw_grid_logX(SDL_Renderer* r,
                           int L, int T, int W, int H, int R, int B,
                           double fmin, double fmax, int nx = 20, int ny = 20) {
    SDL_SetRenderDrawColor(r, 220,220,220,255);
    int x0 = L, x1 = W - R;
    int y0 = T, y1 = H - B;

    // عمودی‌ها در لگاریتم: 20 تقسیم مساوی در log10(f)
    double lmin = std::log10(fmin);
    double lmax = std::log10(fmax);
    for (int i = 1; i < nx; ++i) {
        double lf = lmin + (lmax - lmin) * (double)i / nx;
        int x = x0 + int((lf - lmin) * (x1 - x0) / (lmax - lmin) + 0.5);
        SDL_RenderDrawLine(r, x, y0, x, y1);
    }
    // افقی‌ها خطی (20 قسمت)
    for (int j = 1; j < ny; ++j) {
        int y = y1 - (y1 - y0) * j / ny;
        SDL_RenderDrawLine(r, x0, y, x1, y);
    }
}

static void draw_editor_grid(SDL_Renderer* r, int w, int h, int spacing);

std::vector<double> parse_values_only(const std::string& path) {
    std::setlocale(LC_NUMERIC, "C");

    std::ifstream in(path);
    std::vector<double> vals;
    if (!in.is_open()) return vals;

    std::string line;
    while (std::getline(in, line)) {
        const char* p = line.c_str();
        char* end = nullptr;
        // اولین عدد در خط
        double v = std::strtod(p, &end);
        if (end != p) {
            vals.push_back(v);
        }
    }
    return vals;
}



std::string getElementValueFromUser(const std::string &prompt) {
    char tempPath[MAX_PATH];
    GetTempPathA(MAX_PATH, tempPath);
    // نام فایل یکتا با GetTickCount64 یا استفاده از timestamp
    std::string filePath = std::string(tempPath) + "element_value_" + std::to_string(GetTickCount()) + ".txt";

    std::string scriptPath = std::string(tempPath) + "inputbox.vbs";
    std::ofstream script(scriptPath);
    script << "dim value\n";
    script << "value = InputBox(\"" << prompt << "\", \"value\")\n";
    script << "set fso = CreateObject(\"Scripting.FileSystemObject\")\n";
    script << "set file = fso.CreateTextFile(\"" << filePath << "\", True)\n";
    script << "file.Write value\n";
    script << "file.Close\n";
    script.close();

    // اجرای wscript و صبر برای پایان پردازش (با ShellExecuteEx)
    SHELLEXECUTEINFOA sei = {0};
    sei.cbSize = sizeof(sei);
    sei.fMask = SEE_MASK_NOCLOSEPROCESS;
    sei.lpVerb = "open";
    sei.lpFile = "wscript.exe";
    sei.lpParameters = scriptPath.c_str();
    sei.nShow = SW_HIDE;

    if (ShellExecuteExA(&sei) && sei.hProcess) {
        WaitForSingleObject(sei.hProcess, INFINITE);
        CloseHandle(sei.hProcess);
    } else {
        // اگر ShellExecuteEx شکست خورد، fallback: ساده اجرا کن و کمی تاخیر بده
        ShellExecuteA(NULL, "open", "wscript.exe", scriptPath.c_str(), NULL, SW_HIDE);
        Sleep(300);
    }

    // حالا مطمئنیم فایل ساخته شده — با timeout آن را بخوانیم
    for (int i = 0; i < 200; ++i) {
        std::ifstream result(filePath);
        if (result.is_open()) {
            std::string value;
            std::getline(result, value);
            result.close();
            return value;
        }
        Sleep(10);
    }

    return std::string(); // timeout یا خطا
}

// پنجرهٔ ساده برای انتخاب نوع تحلیل با VBScript MsgBox:
// Yes => Transient, No => AC Sweep, Cancel => Cancel
int chooseAnalysisMode() {
    char tempPath[MAX_PATH];
    GetTempPathA(MAX_PATH, tempPath);
    std::string scriptPath = std::string(tempPath) + "choose_analysis.vbs";

    std::ofstream script(scriptPath);
    script <<
R"(option explicit
dim res
res = InputBox("Select analysis type:" & vbCrLf & "1 = Transient" & vbCrLf & "2 = AC Sweep" & vbCrLf & "3 = Cancel", "Run")
WScript.StdOut.Write res
)";
    script.close();

    SECURITY_ATTRIBUTES sa{sizeof(SECURITY_ATTRIBUTES), nullptr, TRUE};
    HANDLE readPipe, writePipe;
    CreatePipe(&readPipe, &writePipe, &sa, 0);
    SetHandleInformation(readPipe, HANDLE_FLAG_INHERIT, 0);

    STARTUPINFOA si{}; si.cb = sizeof(si);
    si.dwFlags |= STARTF_USESTDHANDLES;
    si.hStdOutput = writePipe;
    si.hStdError  = writePipe;

    PROCESS_INFORMATION pi{};
    std::string cmd = "wscript.exe //Nologo \"" + scriptPath + "\"";
    if (!CreateProcessA(nullptr, (LPSTR)cmd.c_str(), nullptr, nullptr, TRUE, 0, nullptr, nullptr, &si, &pi)) {
        CloseHandle(writePipe); CloseHandle(readPipe);
        return 0;
    }
    CloseHandle(writePipe);

    char buf[16] = {0};
    DWORD readBytes = 0;
    ReadFile(readPipe, buf, sizeof(buf)-1, &readBytes, nullptr);
    CloseHandle(readPipe);
    WaitForSingleObject(pi.hProcess, INFINITE);
    CloseHandle(pi.hThread); CloseHandle(pi.hProcess);

    return atoi(buf); // مستقیماً 1 یا 2 یا 3 برمی‌گرداند
}



bool capture_stdout_to_file(const std::string& path, const std::function<void()>& fn) {
    // 1) ذخیره‌ی stdout فعلی
    int saved_fd = _dup(_fileno(stdout));
    if (saved_fd == -1) return false;

    // 2) بازکردن فایل خروجی
    FILE* new_out = fopen(path.c_str(), "w");
    if (!new_out) {
        _dup2(saved_fd, _fileno(stdout));
        _close(saved_fd);
        return false;
    }

    // 3) هدایت stdout به فایل
    _dup2(_fileno(new_out), _fileno(stdout));

    // 4) اجرای کاری که روی stdout می‌نویسه (phase1_exec)
    fn();
    fflush(stdout);

    // 5) بازگردانی stdout
    _dup2(saved_fd, _fileno(stdout));
    _close(saved_fd);
    fclose(new_out);
    return true;
}


std::vector<std::pair<double,double>> parse_tran_file(const std::string& path) {
    // اطمینان از اعشار با نقطه
    std::setlocale(LC_NUMERIC, "C");

    std::ifstream in(path);
    std::vector<std::pair<double,double>> data;
    if (!in.is_open()) return data;

    std::string line;
    while (std::getline(in, line)) {
        // خط‌هایی مثل "time   V(n)" یا هدر رو نادیده نمی‌گیریم؛
        // به‌جاش همه‌جا دنبال عدد می‌گردیم:
        const char* p = line.c_str();
        char* end = nullptr;
        std::vector<double> nums;
        while (*p) {
            // strtod اگر عدد نبود، end == p برمی‌گرده
            double v = std::strtod(p, &end);
            if (end == p) {
                ++p; // یک کاراکتر جلو برو و دوباره تلاش کن
            } else {
                nums.push_back(v);
                p = end; // از بعدِ عدد ادامه بده
            }
        }
        // اگر حداقل دو عدد در خط بود، اولی = زمان، دومی = مقدار
        if (nums.size() >= 2) {
            data.emplace_back(nums[0], nums[1]);
        }
    }
    return data;
}





using namespace std;
const int SNAP_RADIUS = 5; // snap distance in pixels for wires
int MAGNET_RADIUS = 12;       // جذب به نزدیک‌ترین نود
const int NODE_HIT_RADIUS = 10; // برای کلیک/راست‌کلیک
const int HIT_PAD = 8;          // پدینگ باکس برخورد قطعات
bool showNodeIds = true;

// Union-Find for merging connected nodes
std::vector<int> parent;


static bool snapToNearestNode(int x, int y, int& sx, int& sy) ;
int findParent(int x) {
    if (parent[x] != x) parent[x] = findParent(parent[x]);
    return parent[x];
}

void unionNodes(int a, int b) {
    a = findParent(a);
    b = findParent(b);
    if (a != b) parent[b] = a;
}



void rect(SDL_Renderer* Renderer, int x, int y, int w, int h, int R, int G, int B, int A, int fill_boolian);
void texture(SDL_Renderer* renderer, int x, int y, string address, int w, int h);
void draw_ground(SDL_Renderer* renderer, int x, int y);
void runTransient();
void draw_resistor(SDL_Renderer* renderer, int x, int y, int angle =0);
void draw_capacity(SDL_Renderer* renderer, int x, int y , int angle = 0);
void draw_diode(SDL_Renderer* renderer, int x, int y, int angle = 0);
void draw_voltage(SDL_Renderer* renderer, int x, int , int angle = 0);
void draw_current_source(SDL_Renderer* renderer, int x, int y, int angle = 0);
void draw_inductor(SDL_Renderer* renderer, int x, int y, int angle = 0);
void exportToFile(const std::string& path   ) ;
static void clearAll();                        // NEW
static std::string ask_user_for_save_path();   // NEW

const int GRID_SIZE = 20;
int snap(int value) {
    return (value + GRID_SIZE / 2) / GRID_SIZE * GRID_SIZE;
}
// --- grid-aligned lengths (multiples of GRID_SIZE)
static inline int LEN_R() { return 3 * GRID_SIZE; }  // 60
static inline int LEN_C() { return 2 * GRID_SIZE; }  // 40
static inline int LEN_D() { return 3 * GRID_SIZE; }  // 60
static inline int LEN_V() { return 4 * GRID_SIZE; }  // 80
static inline int LEN_I() { return 4 * GRID_SIZE; }  // 80
static inline int LEN_L() { return 4 * GRID_SIZE; }  // 80
int current_rotation = 0;


enum CursorMode {
    NORMAL,
    PLACING_RESISTOR,
    PLACING_CAPACITOR,
    PLACING_DIODE,
    PLACING_VOLTAGE,
    PLACING_CURRENT_SOURCE,
    PLACING_INDUCTOR,
    PLACING_GROUND,
    PLACING_LABEL_B,
    DRAWING_WIRE,
    DELETE_SELECT   // جدید
};


CursorMode cursorMode = NORMAL;

struct Resistor {
    int x, y;
    int angle = 0;  // زاویه پیش‌فرض 0 درجه
};
vector<Resistor> resistors;

struct Capacitor { int x, y; };
vector<Capacitor> capacitors;

struct Diode { int x, y; };
vector<Diode> diodes;

struct Voltage { int x, y; };
vector<Voltage> voltages;


struct CurrentSrc { int x, y; };
vector<CurrentSrc> currents;


struct Inductor { int x, y; };
vector<Inductor> inductors;

struct Ground { int x, y; };
vector<Ground> grounds;
std::vector<int> ground_node_ids;

struct Wire { int x1, y1, x2, y2; };
vector<Wire> wires;
bool wire_in_progress = false;
int wire_start_x = 0, wire_start_y = 0;

struct LabelB { int x, y, nodeId; };
std::vector<LabelB> labelsB;
int labelB_session_count = 0;  
int activeBGroup = -1;
std::vector<std::vector<int>> labelB_groups; // هر ورودی = گره‌های یک جلسه‌ی B


// --- Node و Element ---
struct Node {
    int id;
    int x, y; // مختصات Grid
};
vector<Node> nodes;



struct Element {
    string type; // "R", "C", "L", "V", "I", "D"
    int node1, node2; // اندیس گره‌ها
    double value;     // مقدار المان
};

// داده‌ی شکل موج برای منابع V/I
// داده‌ی شکل موج برای منابع V/I
struct SourceWF {
    bool isSin   = false;
    bool isPulse = false;
    std::string off, amp, freq;  // برای SIN و PULSE
    std::string pulseType;       // "Step" | "Square" | "Triangular" | "Delta"
};



vector<Element> elements;


// فاصله نقطه تا پاره‌خط (px,py) تا (x1,y1)-(x2,y2)
static double pointToSegmentDist(int px, int py, int x1, int y1, int x2, int y2) {
    double vx = x2 - x1, vy = y2 - y1;
    double wx = px - x1, wy = py - y1;
    double c1 = vx*wx + vy*wy;
    if (c1 <= 0) return std::hypot(px - x1, py - y1);
    double c2 = vx*vx + vy*vy;
    if (c2 <= c1) return std::hypot(px - x2, py - y2);
    double t = c1 / c2;
    double projx = x1 + t * vx, projy = y1 + t * vy;
    return std::hypot(px - (int)projx, py - (int)projy);
}

// برمی‌گرداند index المنت در elements یا -1 اگر چیزی زیر ماوس نبود
static int pickElementAt(int mx, int my, int hitPad = 8) {
    // از آخر به اول تا «آخرین رسم‌شده» اولویت بگیرد
    for (int i = (int)elements.size() - 1; i >= 0; --i) {
        const auto& el = elements[i];
        // گرفتن مختصات دو گرهٔ این المنت
        int x1 = nodes[el.node1].x, y1 = nodes[el.node1].y;
        int x2 = nodes[el.node2].x, y2 = nodes[el.node2].y;
        double d = pointToSegmentDist(mx, my, x1, y1, x2, y2);
        if (d <= hitPad) return i;
    }
    return -1;
}

bool del_select_active = false; // آیا در حال درگ راست برای رسم مستطیل هستیم؟
int sel_x0 = 0, sel_y0 = 0;     // نقطه شروع (راست‌کلیک down)
int sel_x1 = 0, sel_y1 = 0;     // نقطه جاری/پایان


// پیدا کردن یا ساخت گره جدید
int findOrCreateNode(int x, int y) {

    x = snap(x);
    y = snap(y);
    // Snap to existing node if within SNAP_RADIUS
    for (int i = 0; i < (int)nodes.size(); i++) {
        int dx = nodes[i].x - x;
        int dy = nodes[i].y - y;
        if (dx*dx + dy*dy <= SNAP_RADIUS*SNAP_RADIUS) {
            // Found a close enough node, snap to it
            x = nodes[i].x;
            y = nodes[i].y;
            return i;
        }
    }

    for (int i = 0; i < (int)nodes.size(); i++) {
        if (abs(nodes[i].x - x) < GRID_SIZE &&
            abs(nodes[i].y - y) < GRID_SIZE) {
            return i; // گره موجود پیدا شد
        }
    }
    Node newNode { (int)nodes.size(), x, y };
    nodes.push_back(newNode);
    parent.push_back(nodes.size() - 1); // init union-find parent for new node
    return newNode.id;
}

// افزودن المان جدید
void addElement(string type, int x1, int y1, int x2, int y2, double value) {
    int n1 = findOrCreateNode(x1, y1);
    int n2 = findOrCreateNode(x2, y2);
    elements.push_back({type, n1, n2, value});
}



const int SCREEN_WIDTH = 1400;

const int CONNECTOR_SIZE = 5; // اندازه مربع اتصال
const double WIRE_RESISTANCE = 1e-6; // مقاومت خیلی کوچک برای سیم
const int SCREEN_HEIGHT = 760;
//int resistor_offset_x = 0;
// --- wire preview ---
int preview_x2 = 0, preview_y2 = 0;
bool preview_has = false;

void draw_connector(SDL_Renderer* renderer, int cx, int cy) {
    SDL_Rect conn = { cx - CONNECTOR_SIZE/2, cy - CONNECTOR_SIZE/2, CONNECTOR_SIZE, CONNECTOR_SIZE };
    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255); // قرمز
    SDL_RenderFillRect(renderer, &conn);
}

// طول‌های مضرب GRID
static const int RES_LEN = 3 * GRID_SIZE;   // 60 وقتی GRID=20

// گرد کردن به int
static inline int roundi(double v){ return (int)std::round(v); }

// از مختصات سرِ شماره ۱ (روی گرید) و زاویه، مبدأ رسم قبل از دوران را بساز
static inline void originFromAnchorT1(int ax, int ay, int angle, int L, int &ox, int &oy){
    double rad = angle * M_PI / 180.0;
    int cx = ax + roundi((L/2.0) * std::cos(rad));
    int cy = ay + roundi((L/2.0) * std::sin(rad));
    ox = cx - L/2;   // مبدأ local قبل از دوران
    oy = cy;
}

// سرِ دوم از روی سرِ اول و زاویه
static inline void end2FromAnchorT1(int ax, int ay, int angle, int L, int &bx, int &by){
    double rad = angle * M_PI / 180.0;
    bx = ax + roundi(L * std::cos(rad));
    by = ay + roundi(L * std::sin(rad));
}


void show_plot_window(const std::vector<std::pair<double,double>>& data,
                      const std::string& title = "TRAN Plot") {
    if (data.size() < 2) return;

    const int W = 800, H = 480;
    const int L = 60, R = 20, T = 20, B = 40; // مارجین‌های نمودار

    SDL_Window* w = SDL_CreateWindow(title.c_str(),
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, W, H, SDL_WINDOW_SHOWN);
    if (!w) return;
    SDL_Renderer* r = SDL_CreateRenderer(w, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!r) { SDL_DestroyWindow(w); return; }

    // محدودهٔ داده‌ها
    double tmin = data.front().first, tmax = data.back().first;
    for (auto &p : data) { tmin = std::min(tmin,p.first); tmax = std::max(tmax,p.first); }
    double vmin = data.front().second, vmax = data.front().second;
    for (auto &p : data) { vmin = std::min(vmin,p.second); vmax = std::max(vmax,p.second); }
    if (tmax == tmin) tmax = tmin + 1e-12;
    if (vmax == vmin) { vmax = vmin + 1.0; } // جلوگیری از تقسیم بر صفر

    auto X = [&](double t)->int {
        return L + int((t - tmin) * (W - L - R) / (tmax - tmin) + 0.5);
    };
    auto Y = [&](double v)->int {
        return H - B - int((v - vmin) * (H - T - B) / (vmax - vmin) + 0.5);
    };
    CursorUI cu;

    // حلقهٔ کوچک رویداد/رسم (بستن با Esc یا کلید Q یا ضربدر)
    bool quit = false;
    while (!quit) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) quit = true;
            if (e.type == SDL_KEYDOWN) {
                if (e.key.keysym.sym == SDLK_ESCAPE || e.key.keysym.sym == SDLK_q) quit = true;
                if (e.key.keysym.sym == SDLK_c) cu.enabled = !cu.enabled; // کلید میانبر
            }
            if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
                // کلیک روی آیکون → toggle
                if (point_in_rect(e.button.x, e.button.y, cu.icon)) {
                    cu.enabled = !cu.enabled;
                } else if (cu.enabled) {
                    // اگر داخل ناحیهٔ نمودار کلیک شد → نزدیک‌ترین نقطه را انتخاب کن
                    int mx = e.button.x, my = e.button.y;
                    // ناحیهٔ داخل نمودار
                    int plotL=L, plotR=W-R, plotT=T, plotB=H-B;
                    if (mx >= plotL && mx <= plotR && my >= plotT && my <= plotB) {
                        // جستجوی نزدیک‌ترین x
                        // تبدیل پیکسل به زمانِ متناظرِ کلیک
                        double t_click = tmin + (tmax - tmin) * double(mx - (L)) / double((W - R) - (L));
                        // پیدا کردن دو نمونه‌ی قبل/بعد
                        size_t i = bracket_index(data, t_click);
                        double t0 = data[i].first,  y0 = data[i].second;
                        double t1 = data[i+1].first, y1 = data[i+1].second;

                        // میانیابی خطی روی y(t)
                        double alpha = (t1 == t0) ? 0.0 : (t_click - t0) / (t1 - t0);
                        double y_click = y0 + alpha * (y1 - y0);

                        // کراس‌هیر روی نقطه‌ی میانیابی‌شده
                        cu.px = X(t_click);
                        cu.py = Y(y_click);
                        cu.hasPoint = true;

                        // نمایش دقیق در عنوان پنجره
                        char ttl[256];
                        std::snprintf(ttl, sizeof(ttl), "TRAN Cursor: t=%g, y=%g", t_click, y_click);
                        SDL_SetWindowTitle(w, ttl);

                    }
                }
            }
        }


        SDL_SetRenderDrawColor(r, 255,255,255,255); SDL_RenderClear(r);
        draw_grid_linearXY(r, L, T, W, H, R, B, 20, 20);

        // محور‌ها
        SDL_SetRenderDrawColor(r, 0,0,0,255);
        SDL_RenderDrawLine(r, L, H-B, W-R, H-B); // محور t
        SDL_RenderDrawLine(r, L, T,   L,   H-B); // محور v

        // نمودار (Polyline)
        SDL_SetRenderDrawColor(r, 0, 102, 204, 255);
        for (size_t i = 1; i < data.size(); ++i) {
            SDL_RenderDrawLine(r, X(data[i-1].first), Y(data[i-1].second),
                                  X(data[i].first),   Y(data[i].second));
        }
        // نقاط (اختیاری: نمایش نقطه‌های کوچک)
        for (auto &p : data) {
            SDL_Rect dot{ X(p.first)-1, Y(p.second)-1, 3, 3 };
            SDL_RenderFillRect(r, &dot);
        }
        // آیکون حالت نشانگر
        draw_cursor_icon(r, cu.icon, cu.enabled);

        // اگر نقطه‌ای انتخاب شده، کراس‌هیر و نقطه
        if (cu.enabled && cu.hasPoint) {
            SDL_SetRenderDrawColor(r, 200,0,0,255);
            SDL_RenderDrawLine(r, cu.px, T, cu.px, H-B);        // عمودی تا محور x
            SDL_RenderDrawLine(r, L, cu.py, W-R, cu.py);        // افقی تا محور y
            SDL_Rect dot{ cu.px-2, cu.py-2, 5, 5 };
            SDL_RenderFillRect(r, &dot);
        }


        SDL_RenderPresent(r);
    }

    SDL_DestroyRenderer(r);
    SDL_DestroyWindow(w);
}

// ایندکس عنصر در vector<elements> → مشخصات سینوسی
static std::map<int, SourceWF> srcWFByIdx;

// کمکی: تشخیص و پارس "SIN(...)" از یک رشته‌ی ورودی
static bool parseSIN(const std::string& s,
                     std::string& off, std::string& amp, std::string& freq) {
    static const std::regex re(
        R"(^\s*(?:SIN|sin)\s*\(\s*([^\s\)]+)\s+([^\s\)]+)\s+([^\s\)]+)\s*\)\s*$)",
        std::regex::icase
    );
    std::smatch m;
    if (std::regex_match(s, m, re) && m.size() == 4) {
        off  = m[1], amp = m[2], freq = m[3];
        return true;
    }
    return false;
}
// تشخیص "PULSE <Type> (<off> <amp> <freq>)"  یا  "PULSE Delta"
static bool parsePULSE(const std::string& s,
                       std::string& type, std::string& off, std::string& amp, std::string& freq) {
    static const std::regex reFull(
        R"(^\s*(?:PULSE|pulse)\s+([A-Za-z]+)\s*\(\s*([^\s\)]+)\s+([^\s\)]+)\s+([^\s\)]+)\s*\)\s*$)",
        std::regex::icase
    );
    std::smatch m;
    if (std::regex_match(s, m, reFull) && m.size()==5) {
        type = m[1]; off = m[2]; amp = m[3]; freq = m[4];
        return true;
    }
    static const std::regex reDelta(
        R"(^\s*(?:PULSE|pulse)\s+(Delta)\s*$)",
        std::regex::icase
    );
    if (std::regex_match(s, m, reDelta)) {
        type = "Delta"; off = "0"; amp = "0"; freq = "0";
        return true;
    }
    return false;
}


// پارس فایل AC: انتظار هر خط حداقل دو عدد: freq  value
static std::vector<std::pair<double,double>> parse_ac_pairs(const std::string& path) {
    std::setlocale(LC_NUMERIC, "C");
    std::ifstream in(path);
    std::vector<std::pair<double,double>> data;
    if (!in.is_open()) return data;

    std::string line;
    while (std::getline(in, line)) {
        const char* p = line.c_str();
        char* end = nullptr;
        std::vector<double> nums;
        while (*p) {
            double v = std::strtod(p, &end);
            if (end == p) ++p; else { nums.push_back(v); p = end; }
        }
        if (nums.size() >= 2) data.emplace_back(nums[0], nums[1]);
    }
    return data;
}

// نمایش نمودار Bode Magnitude (X=log10(f))
static void show_bode_window(const std::vector<std::pair<double,double>>& data,
                             bool use_dB, const std::string& title = "AC Plot") {
    if (data.size() < 2) return;

    const int W=900, H=520;
    const int L=60, R=20, T=20, B=40;

    SDL_Window* w = SDL_CreateWindow(title.c_str(),
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, W, H, SDL_WINDOW_SHOWN);
    if (!w) return;
    SDL_Renderer* r = SDL_CreateRenderer(w, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!r) { SDL_DestroyWindow(w); return; }

    // محدودهٔ فرکانس و مقدار
    double fmin = data.front().first, fmax = data.front().first;
    for (auto &p : data) { fmin = std::min(fmin, p.first); fmax = std::max(fmax, p.first); }
    if (fmin <= 0) fmin = 1e-3;
    if (fmax <= fmin) fmax = fmin * 10.0;

    auto toMag = [&](double m){ return use_dB ? 20.0*std::log10(std::max(m,1e-300)) : m; };
    double ymin = toMag(data.front().second), ymax = ymin;
    for (auto &p : data) { ymin = std::min(ymin, toMag(p.second)); ymax = std::max(ymax, toMag(p.second)); }
    if (ymax == ymin) ymax = ymin + 1.0;

    auto X = [&](double f)->int {
        double lf = std::log10(f), lmin = std::log10(fmin), lmax = std::log10(fmax);
        return L + int((lf - lmin) * (W - L - R) / (lmax - lmin) + 0.5);
    };
    auto Y = [&](double v)->int {
        double vv = toMag(v);
        return H - B - int((vv - ymin) * (H - T - B) / (ymax - ymin) + 0.5);
    };
    CursorUI cu;


    bool quit=false;
    while (!quit) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type==SDL_QUIT) quit=true;
            if (e.type==SDL_KEYDOWN && (e.key.keysym.sym==SDLK_ESCAPE || e.key.keysym.sym==SDLK_q)) quit=true;
            if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
                if (point_in_rect(e.button.x, e.button.y, cu.icon)) {
                    cu.enabled = !cu.enabled;
                } else if (cu.enabled) {
                    int mx = e.button.x, my = e.button.y;
                    int plotL=L, plotR=W-R, plotT=T, plotB=H-B;
                    if (mx >= plotL && mx <= plotR && my >= plotT && my <= plotB) {
                        int mx = e.button.x, my = e.button.y;
                        int plotL=L, plotR=W-R, plotT=T, plotB=H-B;
                        if (mx >= plotL && mx <= plotR && my >= plotT && my <= plotB) {
                            // پیکسل → فرکانس (محور لگاریتمی)
                            double lmin = std::log10(fmin), lmax = std::log10(fmax);
                            double lf_click = lmin + (lmax - lmin) * double(mx - plotL) / double(plotR - plotL);
                            double f_click = std::pow(10.0, lf_click);

                            // اندیس قبل/بعد نسبت به فرکانس
                            size_t i = bracket_index(data, f_click);
                            double f0 = data[i].first,  m0 = data[i].second; // magnitude خطی در فایل
                            double f1 = data[i+1].first, m1 = data[i+1].second;

                            // میانیابی خطی روی مقدار «خطی»، نه dB
                            double alpha = (f1 == f0) ? 0.0 : (f_click - f0) / (f1 - f0);
                            double mag_lin = m0 + alpha * (m1 - m0);

                            // تبدیل برای نمایش (اگر dB انتخاب شده)
                            double mag_disp = use_dB ? 20.0 * std::log10(std::max(mag_lin, 1e-300)) : mag_lin;

                            // پیکسل نقطهٔ میانیابی‌شده
                            cu.px = X(f_click);
                            cu.py = Y(mag_lin);  // توجه: Y() خودش با dB/lin هماهنگ شده چون از toMag استفاده کردی
                            cu.hasPoint = true;

                            // عنوان پنجره: هر دو حالت lin و dB را نشان بده
                            char ttl[256];
                            if (use_dB) {
                                std::snprintf(ttl, sizeof(ttl), "AC Cursor: f=%g Hz, mag=%g (lin), %g dB",
                                              f_click, mag_lin, mag_disp);
                            } else {
                                std::snprintf(ttl, sizeof(ttl), "AC Cursor: f=%g Hz, mag=%g", f_click, mag_lin);
                            }
                            SDL_SetWindowTitle(w, ttl);
                        }

                    }
                }
            }

        }
        // پس‌زمینه
        SDL_SetRenderDrawColor(r, 255,255,255,255); SDL_RenderClear(r);
        draw_grid_logX(r, L, T, W, H, R, B, fmin, fmax, 20, 20);

        // محور‌ها
        SDL_SetRenderDrawColor(r, 0,0,0,255);
        SDL_RenderDrawLine(r, L, H-B, W-R, H-B); // فرکانس
        SDL_RenderDrawLine(r, L, T,   L,   H-B); // مقدار

        // نمودار
        SDL_SetRenderDrawColor(r, 0, 102, 204, 255);
        for (size_t i=1;i<data.size();++i) {
            SDL_RenderDrawLine(r, X(data[i-1].first), Y(data[i-1].second),
                                  X(data[i].first),   Y(data[i].second));
        }
        for (auto &p: data) {
            SDL_Rect dot{ X(p.first)-1, Y(p.second)-1, 3, 3 };
            SDL_RenderFillRect(r, &dot);
        }
        draw_cursor_icon(r, cu.icon, cu.enabled);
        if (cu.enabled && cu.hasPoint) {
            SDL_SetRenderDrawColor(r, 200,0,0,255);
            SDL_RenderDrawLine(r, cu.px, T, cu.px, H-B);
            SDL_RenderDrawLine(r, L, cu.py, W-R, cu.py);
            SDL_Rect dot{ cu.px-2, cu.py-2, 5, 5 };
            SDL_RenderFillRect(r, &dot);
        }

        SDL_RenderPresent(r);
    }
    SDL_DestroyRenderer(r); SDL_DestroyWindow(w);
}



void runTransient() {
    try {
        // 1) Export the current circuit to temp file
        exportToFile("temp_run.txt");

        // 2) Load it into Phase 1 (phase1_load_file clears prior state)
        phase1_load_file("temp_run.txt");

        // 3) Ask user for transient parameters
        std::string params = getElementValueFromUser(
            "Enter TRAN params: Step Stop Start  (like: 0.01 0.1 0)"
        );
        if (params.empty()) {
            MessageBoxA(NULL, "No parameters entered.", "Run aborted", MB_OK | MB_ICONWARNING);
            return;
        }

        std::istringstream iss(params);
        double step = 0.0, stop = 0.0, start = 0.0;
        if (!(iss >> step >> stop >> start)) {
            MessageBoxA(NULL,
                        "Invalid format. Use: <time_step> <total_time> <start_time>",
                        "Error", MB_OK | MB_ICONERROR);
            return;
        }

        // خروجی‌ها: V(node) یا I(component)
        std::string outputs = getElementValueFromUser(
            "Outputs: V(node_number) or I(component_name) (default = V(node))"
        );
        if (outputs.empty()) {
            MessageBoxA(NULL, "No outputs entered. Example: V(2) or I(V1)",
                        "Run aborted", MB_OK | MB_ICONWARNING);
            return;
        }

        // 5) Build and execute the .print TRAN command
        // 5) Build the .print TRAN command
        std::ostringstream cmd;
        cmd << ".print TRAN " << step << " " << stop << " " << start << " 1 " << outputs;

        // 6) اجرا و هدایت خروجی کنسول به فایل، سپس پارس و رسم
        std::string outPath = std::string("tran_out.txt");
        bool ok = capture_stdout_to_file(outPath, [&](){
            phase1_exec(cmd.str());   // همان کال فعلی شما
        });

        if (!ok) {
            MessageBoxA(NULL, "Failed to capture output.", "Plot", MB_OK | MB_ICONERROR);
            return;
        }

        // 7) خواندن فایل و رسم نمودار
        // 7) خواندن فایل و ساخت داده‌های (t, y)
        auto values = parse_values_only(outPath);
        if (values.size() < 2) {
            MessageBoxA(NULL, "No time-series data parsed.\nFile has only one or zero values.", "Plot", MB_OK | MB_ICONWARNING);
            return;
        }

        // تعداد نمونهٔ انتظار رفته از روی start/stop/step
        size_t n_expected = 0;
        if (step > 0 && stop >= start) {
            n_expected = static_cast<size_t>(std::floor((stop - start) / step + 1e-12)) + 1;
        }
        size_t n = values.size();
        if (n_expected > 0) n = std::min(n, n_expected);

        // ساخت داده‌های (t, y)
        std::vector<std::pair<double,double>> data;
        data.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            double t = start + i * step;
            data.emplace_back(t, values[i]);
        }

        show_plot_window(data, "TRAN Plot: " + outputs);


    } catch (const std::exception& ex) {
        MessageBoxA(NULL, ex.what(), "Run error", MB_OK | MB_ICONERROR);
        return;
    }
}

void runACSweep() {
    try {
        // 1) Export و Load
        exportToFile("temp_run.txt");
        phase1_load_file("temp_run.txt");

        // 2) پارامترهای AC: نوع، نقاط، fstart، fstop
        std::string stype = getElementValueFromUser("AC sweep type (dec or lin):");
        if (stype.empty()) return;
        for (auto &ch : stype) ch = (char)tolower((unsigned char)ch);
        if (stype != "dec" && stype != "lin") {
            MessageBoxA(NULL, "Type must be 'dec' or 'lin'", "AC", MB_OK | MB_ICONERROR);
            return;
        }

        std::string spts  = getElementValueFromUser(stype=="dec" ?
            "Points per decade (e.g. 20):" : "Total points (e.g. 200):");
        if (spts.empty()) return;
        int points = atoi(spts.c_str());
        if (points < 2) { MessageBoxA(NULL, "Points must be >= 2", "AC", MB_OK|MB_ICONERROR); return; }

        std::string sfstart = getElementValueFromUser("Start frequency (Hz), e.g. 10:");
        if (sfstart.empty()) return;
        double fstart = atof(sfstart.c_str());
        std::string sfstop  = getElementValueFromUser("Stop frequency (Hz), e.g. 1e6:");
        if (sfstop.empty()) return;
        double fstop  = atof(sfstop.c_str());
        if (fstart <= 0 || fstop <= fstart) {
            MessageBoxA(NULL, "Invalid frequency range", "AC", MB_OK|MB_ICONERROR); return;
        }

        // 3) خروجی مورد نظر
        std::string outputs = getElementValueFromUser("Output(s): e.g. V(2) or I(V1)");
        if (outputs.empty()) {
            MessageBoxA(NULL, "No outputs entered.", "AC", MB_OK | MB_ICONWARNING);
            return;
        }

        // 4) دسی‌بل یا خطی؟
        std::string ymode = getElementValueFromUser("Magnitude mode: 'db' or 'lin' (default lin):");
        bool use_dB = false;
        if (!ymode.empty()) {
            for (auto &c: ymode) c = (char)tolower((unsigned char)c);
            use_dB = (ymode == "db");
        }

        // 5) ساخت دستور و اجرا
        std::ostringstream cmd;
        // فرمت با کلاس Analysis فاز1 هماهنگ است: .print AC <dec|lin> <points> <fstart> <fstop> <outputs...>
        cmd << ".print AC " << stype << " " << points << " " << fstart << " " << fstop << " " << outputs;

        std::string outPath = "ac_out.txt";
        bool ok = capture_stdout_to_file(outPath, [&](){
            phase1_exec(cmd.str());
        });
        if (!ok) { MessageBoxA(NULL, "Failed to capture output.", "AC", MB_OK | MB_ICONERROR); return; }

        // 6) پارس و رسم
        auto data = parse_ac_pairs(outPath);
        if (data.size() < 2) {
            MessageBoxA(NULL, "No AC data parsed. Check .print format.", "AC", MB_OK | MB_ICONWARNING);
            return;
        }
        std::string ttl = std::string("AC Plot (") + (use_dB ? "dB" : "lin") + "): " + outputs;
        show_bode_window(data, use_dB, ttl);

    } catch (const std::exception& ex) {
        MessageBoxA(NULL, ex.what(), "AC error", MB_OK | MB_ICONERROR);
        return;
    }
}


static inline bool insideRect(int x, int y, int rx, int ry, int rw, int rh) {
    return x >= rx && x <= rx + rw && y >= ry && y <= ry + rh;
}
static inline void normRect(int x0, int y0, int x1, int y1, int &rx, int &ry, int &rw, int &rh) {
    rx = std::min(x0, x1); ry = std::min(y0, y1);
    rw = std::abs(x1 - x0); rh = std::abs(y1 - y0);
}

// پیدا کردن گره موجود در مختصات (بدون ساخت گره جدید)
int getExistingNodeIdAt(int x, int y) {
    for (int i = 0; i < (int)nodes.size(); ++i) {
        int dx = nodes[i].x - x;
        int dy = nodes[i].y - y;
        if (dx*dx + dy*dy <= SNAP_RADIUS*SNAP_RADIUS) return i;
    }
    return -1;
}


int main(int argc, char* argv[]) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        cerr << "SDL could not initialize! SDL_Error: " << SDL_GetError() << endl;
        return -1;
    }
    if (TTF_Init() == -1) {
        cerr << "TTF could not initialize! TTF_Error: " << TTF_GetError() << endl;
        SDL_Quit();
        return -1;
    }

    SDL_Window* window = SDL_CreateWindow("SDL Text Editor",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);

    if (!window) {
        cerr << "Window could not be created! SDL_Error: " << SDL_GetError() << endl;
        TTF_Quit(); SDL_Quit(); return -1;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!renderer) {
        cerr << "Renderer could not be created! SDL_Error: " << SDL_GetError() << endl;
        SDL_DestroyWindow(window); TTF_Quit(); SDL_Quit(); return -1;
    }

    SDL_Cursor* cursor = SDL_CreateSystemCursor(SDL_SYSTEM_CURSOR_CROSSHAIR);
    SDL_SetCursor(cursor);

    bool quit = false;
    SDL_Event e;

    while (!quit) {
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT)
                quit = true;

            if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
                int mx = e.button.x;
                int my = e.button.y;
                if (mx >= 630 && mx <= 656 && my >= 20 && my <= 48) {
                    if (cursorMode == DELETE_SELECT) {
                        cursorMode = NORMAL;
                        del_select_active = false; // اگر مستطیل نیمه‌کاره بود
                    } else {
                        cursorMode = DELETE_SELECT;
                        del_select_active = false;
                    }
                    // مهم: اینجا return نکن؛ اجازه بده باقی ifها بررسی بشن فقط اگه لازم داری.
                }

                if (cursorMode == DELETE_SELECT) {
                    del_select_active = true;
                    sel_x0 = sel_x1 = snap(e.button.x);
                    sel_y0 = sel_y1 = snap(e.button.y);
                    // توجه: در این حالت، دیگر شرط قرار دادن قطعه/سیم اجرا نشود.
                }

                // NEW BUTTON: (20,20)-(46,48)
                if (cursorMode == NORMAL &&
                    e.button.x >= 20 && e.button.x <= (20+26) &&
                    e.button.y >= 20 && e.button.y <= (20+28)) {
                        clearAll();
                        cursorMode = NORMAL;
                        SDL_ShowCursor(SDL_ENABLE);
                        SDL_SetCursor(cursor);
                    }

                                // SAVE BUTTON: (55,20)-(81,48)
                                if (cursorMode == NORMAL &&
                                    e.button.x >= 55 && e.button.x <= (55+26) &&
                                    e.button.y >= 20 && e.button.y <= (20+28)) {
                                        std::string path = ask_user_for_save_path();
                                        if (!path.empty()) {
                                                exportToFile(path);  // همان خروجی متن فاز۱
                                                MessageBoxA(NULL, ("Saved to:\n" + path).c_str(), "Save", MB_OK | MB_ICONINFORMATION);
                                            }
                                    }

                if (cursorMode == NORMAL && mx >= 700 && mx <= 726 && my >= 20 && my <= 48) {
                    cursorMode = PLACING_RESISTOR;
                    SDL_ShowCursor(SDL_DISABLE);
                }

                else if (cursorMode == NORMAL && mx >= 945 && mx <= 971 && my >= 20 && my <= 48) {
                    cursorMode = PLACING_LABEL_B;
                    labelB_session_count = 0;
                    activeBGroup = (int)labelB_groups.size();
                    labelB_groups.push_back({}); // گروه جدید
                    SDL_ShowCursor(SDL_DISABLE);
                }
                else if (cursorMode == PLACING_LABEL_B && mx >= 20 && mx <= 1100 && my >= 20 && my <= 48) {
                    if (labelB_session_count < 2) {
                        MessageBoxA(NULL, "2 Points needed at least", "B label", MB_OK | MB_ICONWARNING);
                    } else {
                        cursorMode = NORMAL;
                        SDL_ShowCursor(SDL_ENABLE);
                        SDL_SetCursor(cursor);
                        activeBGroup = -1; // گروه فعلی نهایی شد
                    }
                }
                else if (cursorMode == PLACING_RESISTOR) {
                    int cx   = snap(mx);
                    int cy   = snap(my);
                    int left = snap(cx - LEN_R()/2);
                    resistors.push_back({ left, cy, current_rotation });

                    int n1 = findOrCreateNode(left, cy);
                    int n2 = findOrCreateNode(left + LEN_R(), cy);

                    std::string valStr = getElementValueFromUser("enter R value:");
                    double val = atof(valStr.c_str());
                    elements.push_back({ "R", n1, n2, val });

                    cursorMode = NORMAL;
                    SDL_ShowCursor(SDL_ENABLE);
                    SDL_SetCursor(cursor);
                }
                else if (cursorMode == PLACING_LABEL_B) {
                    int bx = snap(mx);
                    int by = snap(my);
                    int nid = findOrCreateNode(bx, by);
                    labelsB.push_back({bx, by, nid});
                    if (activeBGroup >= 0) labelB_groups[activeBGroup].push_back(nid);
                    labelB_session_count++;
                }
                else if (cursorMode == PLACING_CURRENT_SOURCE) {
                    int cx = snap(mx), cy = snap(my);
                    int left = cx - LEN_I()/2;
                    currents.push_back({ left, cy });

                    int n1 = findOrCreateNode(left, cy);
                    int n2 = findOrCreateNode(left + LEN_I(), cy);

                    std::string iStr = getElementValueFromUser("enter I: (e.g. 2 or SIN(0 5 1k) or PULSE Square (0 5 1k))");
                    int idx = (int)elements.size();

                    std::string off, amp, freq, ptype;
                    if (parseSIN(iStr, off, amp, freq)) {
                        elements.push_back({ "I", n1, n2, 0.0 });
                        srcWFByIdx[idx] = { true, false, off, amp, freq, "" };
                    } else if (parsePULSE(iStr, ptype, off, amp, freq)) {
                        elements.push_back({ "I", n1, n2, 0.0 });
                        SourceWF wf; wf.isSin=false; wf.isPulse=true; wf.off=off; wf.amp=amp; wf.freq=freq; wf.pulseType=ptype;
                        srcWFByIdx[idx] = wf;
                    } else {
                        double val = atof(iStr.c_str());
                        elements.push_back({ "I", n1, n2, val });
                    }

                    cursorMode = NORMAL; SDL_ShowCursor(SDL_ENABLE); SDL_SetCursor(cursor);
                }

                else if (cursorMode == NORMAL && mx >= 735 && mx <= 761 && my >= 20 && my <= 48) {
                    cursorMode = PLACING_CAPACITOR;
                    SDL_ShowCursor(SDL_DISABLE);
                }
                else if (cursorMode == PLACING_CAPACITOR) {
                    int cx = snap(mx), cy = snap(my);
                    int left = cx - LEN_C()/2;
                    capacitors.push_back({ left, cy });

                    int n1 = findOrCreateNode(left, cy);
                    int n2 = findOrCreateNode(left + LEN_C(), cy);

                    std::string valStr = getElementValueFromUser("enter C value:");
                    double val = atof(valStr.c_str());
                    elements.push_back({ "C", n1, n2, val });

                    cursorMode = NORMAL; SDL_ShowCursor(SDL_ENABLE); SDL_SetCursor(cursor);
                }
                else if (cursorMode == NORMAL && mx >= 770 && mx <= 796 && my >= 20 && my <= 48) {
                    cursorMode = PLACING_DIODE;
                    SDL_ShowCursor(SDL_DISABLE);
                }
                else if (cursorMode == PLACING_DIODE) {
                    int cx = snap(mx), cy = snap(my);
                    int left = cx - LEN_D()/2;
                    diodes.push_back({ left, cy });

                    int n1 = findOrCreateNode(left, cy);
                    int n2 = findOrCreateNode(left + LEN_D(), cy);
                    elements.push_back({ "D", n1, n2, 0 });

                    cursorMode = NORMAL; SDL_ShowCursor(SDL_ENABLE); SDL_SetCursor(cursor);
                }

                else if (cursorMode == NORMAL && mx >= 805 && mx <= 831 && my >= 20 && my <= 48) {
                    cursorMode = PLACING_VOLTAGE;
                    SDL_ShowCursor(SDL_DISABLE);
                }
                else if (cursorMode == PLACING_VOLTAGE) {
                    int cx = snap(mx), cy = snap(my);
                    int left = cx - LEN_V()/2;
                    voltages.push_back({ left, cy });

                    int n1 = findOrCreateNode(left, cy);
                    int n2 = findOrCreateNode(left + LEN_V(), cy);

                    std::string vStr = getElementValueFromUser("enter V: (e.g. 5 or SIN(0 5 1k) or PULSE Square (0 5 1k))");
                    int idx = (int)elements.size();

                    std::string off, amp, freq, ptype;
                    if (parseSIN(vStr, off, amp, freq)) {
                        elements.push_back({ "V", n1, n2, 0.0 });
                        srcWFByIdx[idx] = { true, false, off, amp, freq, "" };
                    } else if (parsePULSE(vStr, ptype, off, amp, freq)) {
                        elements.push_back({ "V", n1, n2, 0.0 });
                        SourceWF wf; wf.isSin=false; wf.isPulse=true; wf.off=off; wf.amp=amp; wf.freq=freq; wf.pulseType=ptype;
                        srcWFByIdx[idx] = wf;
                    } else {
                        double val = atof(vStr.c_str());
                        elements.push_back({ "V", n1, n2, val });
                    }

                    cursorMode = NORMAL; SDL_ShowCursor(SDL_ENABLE); SDL_SetCursor(cursor);
                }


                else if (cursorMode == NORMAL && mx >= 840 && mx <= 866 && my >= 20 && my <= 48) {
                    cursorMode = PLACING_INDUCTOR;
                    SDL_ShowCursor(SDL_DISABLE);
                }
                else if (cursorMode == PLACING_INDUCTOR) {
                    int cx = snap(mx), cy = snap(my);
                    int left = cx - LEN_L()/2;
                    inductors.push_back({ left, cy });

                    int n1 = findOrCreateNode(left, cy);
                    int n2 = findOrCreateNode(left + LEN_L(), cy);

                    std::string valStr = getElementValueFromUser("enter L value:");
                    double val = atof(valStr.c_str());
                    elements.push_back({ "L", n1, n2, val });

                    cursorMode = NORMAL; SDL_ShowCursor(SDL_ENABLE); SDL_SetCursor(cursor);
                }

                else if (cursorMode == NORMAL && mx >= 665 && mx <= 690 && my >= 20 && my <= 48) {
                    cursorMode = DRAWING_WIRE;
                }
                else if (cursorMode == DRAWING_WIRE && mx >= 90 && mx <= 690 && my >= 20 && my <= 48) {
                    cursorMode = NORMAL;
                    SDL_ShowCursor(SDL_ENABLE);
                    SDL_SetCursor(cursor);
                }

                else if (cursorMode == NORMAL && mx >= 875 && mx <= 903 && my >= 20 && my <= 48) {
                    cursorMode = PLACING_GROUND;
                    SDL_ShowCursor(SDL_DISABLE);
                }
                else if (cursorMode == NORMAL && mx >= 910 && mx <= 936 && my >= 20 && my <= 48) {
                    cursorMode = PLACING_CURRENT_SOURCE;
                    SDL_ShowCursor(SDL_DISABLE);
                }
                else if (cursorMode == PLACING_GROUND) {
                    int cx = snap(mx), cy = snap(my);   // مرکز آیکون روی گرید
                    int left = cx - 14;                 // چون آیکون 28px است
                    grounds.push_back({ left, cy });

                    int nodeId = findOrCreateNode(cx, cy); // کانکتور بالای آیکون = مرکز
                    ground_node_ids.push_back(nodeId);

                    cursorMode = NORMAL; SDL_ShowCursor(SDL_ENABLE); SDL_SetCursor(cursor);
                }


                else if (cursorMode == DRAWING_WIRE) {
                    if (!wire_in_progress) {
                        int sx = snap(mx), sy = snap(my);
                        if (snapToNearestNode(sx, sy, sx, sy)) { /* sx,sy به نود چسبید */ }
                        wire_start_x = sx; wire_start_y = sy;

                        wire_in_progress = true;
                    } else {
                        int x2 = snap(mx), y2 = snap(my);
                        int sx=x2, sy=y2;
                        if (snapToNearestNode(x2, y2, sx, sy)) { x2 = sx; y2 = sy; }
                        wires.push_back({ wire_start_x, wire_start_y, x2, y2 });
                        wire_in_progress = false;
                        preview_has = false;

                        {
                            int n1 = findOrCreateNode(wire_start_x, wire_start_y);
                            int n2 = findOrCreateNode(x2, y2);
                            elements.push_back({"R", n1, n2, WIRE_RESISTANCE});
                        }
                    }
                }
                // ناحیه‌ی آیکون Run (هم‌ اندازه‌ی PNG: 26x28 در مختصات 90,20)
                // ناحیه آیکون Run (90,20 با اندازه 26x28)
                // آیکون Run: اول انتخاب نوع تحلیل
                // آیکون Run: انتخاب نوع تحلیل با ورودی عددی
                if (cursorMode == NORMAL &&
                    e.button.x >= 90 && e.button.x <= (90+26) &&
                    e.button.y >= 20 && e.button.y <= (20+28)) {

                    int res = chooseAnalysisMode(); // 1=Transient, 2=AC Sweep, 3=Cancel

                    if (res == 1) {
                        runTransient();
                    } else if (res == 2) {
                        runACSweep();
                    } else {
                        // Cancel یا ورودی اشتباه -> هیچ کاری انجام نمی‌ده
                    }
                    }


            }
            else if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_RIGHT) {
                int mx = e.button.x;
                int my = e.button.y;

                // اول: المنت زیر ماوس؟
                int eid = pickElementAt(mx, my, /*hitPad=*/ std::max(8, GRID_SIZE/2));
                if (eid >= 0) {
                    const auto& el = elements[eid];
                    std::ostringstream ss;
                    ss << elName((size_t)eid)              // مثل R3
                       << "  type=" << el.type             // نوع مثل "R"
                       << "  value=" << el.value;          // مقدار (برای V/I DC هم نمایش می‌دهد)

                    MessageBoxA(NULL, ss.str().c_str(), "Element", MB_OK | MB_ICONINFORMATION);
                    // اگر موج V/I غیر DC هم داری و می‌خواهی نشان بدهی، از map مربوط به موج‌ها بخوانی (srcWFByIdx)
                    // (ایندکس همان eid است که قبلاً ذخیره کرده‌ای.)
                    break;
                }

                // دوم: اگر المنت نبود، می‌تونی نود نزدیک را هم چک کنی (اختیاری):
                // از تابع آمادهٔ خودت استفاده کن:
                int nid = getExistingNodeIdAt(mx, my);  // همین الان در کدت موجود است
                if (nid >= 0) {
                    std::string msg = "Node (" + std::to_string(nodes[nid].x) + "," + std::to_string(nodes[nid].y) + ")";
                    MessageBoxA(NULL, msg.c_str(), "Node", MB_OK | MB_ICONINFORMATION);
                    break;
                }
            }

            else if (e.type == SDL_KEYDOWN) {
                if (e.key.keysym.sym == SDLK_n) showNodeIds = !showNodeIds;
                if (e.key.keysym.sym == SDLK_e) {
    exportToFile("circuit.txt");
                    // یک پیام ساده هم نمایش بده (اختیاری):
                    MessageBoxA(NULL, "Circuit exported to circuit.txt", "Export", MB_OK | MB_ICONINFORMATION);
                }

                const Uint8* state = SDL_GetKeyboardState(nullptr);

                if (state[SDL_SCANCODE_LCTRL] || state[SDL_SCANCODE_RCTRL]) {
                    if (e.key.keysym.sym == SDLK_r) {
                        if (cursorMode == PLACING_RESISTOR) {
                            current_rotation = (current_rotation + 90) % 360;
                           // resistor_offset_x -= 5;
                        }
                    }
                }
            }
            if (cursorMode == DELETE_SELECT && del_select_active && (e.motion.state & SDL_BUTTON_LMASK)) {
                sel_x1 = snap(e.motion.x);
                sel_y1 = snap(e.motion.y);
            }
            else if (e.type == SDL_MOUSEMOTION) {
                if (cursorMode == DRAWING_WIRE && wire_in_progress) {
                    // مختصات ماوس → اسنپ روی گرید
                    int mx = snap(e.motion.x);
                    int my = snap(e.motion.y);

                    preview_x2 = mx;
                    preview_y2 = my;
                    preview_has = true;
                }
            }

            if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_LEFT) {
    if (cursorMode == DELETE_SELECT && del_select_active) {
        del_select_active = false;

        int rx, ry, rw, rh;
        normRect(sel_x0, sel_y0, snap(e.button.x), snap(e.button.y), rx, ry, rw, rh);

        // --- حذف اشیا فقط اگر «کامل» داخل مستطیل‌اند ---

        // 1) wires: هر دو سر داخل => حذف
        {
            std::vector<Wire> keep;
            for (auto &w : wires) {
                bool in1 = insideRect(w.x1, w.y1, rx, ry, rw, rh);
                bool in2 = insideRect(w.x2, w.y2, rx, ry, rw, rh);
                if (!(in1 && in2)) keep.push_back(w);
            }
            wires.swap(keep);
        }

        // 2) remove_element_by_2nodes helper (مثل نسخه‌ای که قبلاً دادم/گذاشتی)
        auto remove_element_by_2nodes = [&](const std::string &type, int nx1, int ny1, int nx2, int ny2) {
            int id1 = getExistingNodeIdAt(nx1, ny1);
            int id2 = getExistingNodeIdAt(nx2, ny2);
            if (id1 < 0 || id2 < 0) return;
            for (size_t i = 0; i < elements.size(); ) {
                bool sameType = (elements[i].type == type);
                bool sameNodes = ((elements[i].node1 == id1 && elements[i].node2 == id2) ||
                                  (elements[i].node1 == id2 && elements[i].node2 == id1));
                if (sameType && sameNodes) elements.erase(elements.begin()+i);
                else ++i;
            }
        };

        // 3) قطعات نمایشی + elements متناظر (مرزهای کامل داخل مستطیل)
        // مقاومت: دو سر نمایشی بعد از چرخش را چک کن (تابع rotated_res_ends مثل قبل)
        // قبلی که 70 ثابت داشت را با این جایگزین کن:
        auto rotated_res_ends = [&](int x, int y, int angle, int &ex1, int &ey1, int &ex2, int &ey2) {
            int L = LEN_R();                  // ← طول واقعی مقاومت (مثلاً 60)
            int px0 = x,      py0 = y;
            int px7 = x + L,  py7 = y;
            int cx  = x + L/2, cy = y;
            double rad = angle * M_PI / 180.0;
            auto rot = [&](int &xx, int &yy){
                int tx = xx - cx, ty = yy - cy;
                int rx = int(tx * cos(rad) - ty * sin(rad));
                int ry = int(tx * sin(rad) + ty * cos(rad));
                xx = cx + rx; yy = cy + ry;
            };
            rot(px0, py0); rot(px7, py7);
            ex1 = px0; ey1 = py0; ex2 = px7; ey2 = py7;
        };


        { // R
            std::vector<Resistor> keep;
            for (auto &r : resistors) {
                int ex1, ey1, ex2, ey2;
                rotated_res_ends(r.x, r.y, r.angle, ex1, ey1, ex2, ey2);
                bool in1 = insideRect(ex1, ey1, rx, ry, rw, rh);
                bool in2 = insideRect(ex2, ey2, rx, ry, rw, rh);
                if (in1 && in2) {
                    // به remove_element_by_2nodes مختصات همان سرها را بده
                    remove_element_by_2nodes("R", ex1, ey1, ex2, ey2);
                } else {
                    keep.push_back(r);
                }
            }
            resistors.swap(keep);
        }
        { // C: از x تا x+LEN_C()
            std::vector<Capacitor> keep;
            for (auto &c : capacitors) {
                int ex1=c.x, ey1=c.y, ex2=c.x + LEN_C(), ey2=c.y;   // ← به‌جای 35
                bool in1=insideRect(ex1,ey1,rx,ry,rw,rh), in2=insideRect(ex2,ey2,rx,ry,rw,rh);
                if (in1 && in2) remove_element_by_2nodes("C", ex1, ey1, ex2, ey2);
                else keep.push_back(c);
            }
            capacitors.swap(keep);
        }
        { // D: از x تا x+LEN_D()
            std::vector<Diode> keep;
            for (auto &d : diodes) {
                int ex1=d.x, ey1=d.y, ex2=d.x + LEN_D(), ey2=d.y;   // ← به‌جای 45
                bool in1=insideRect(ex1,ey1,rx,ry,rw,rh), in2=insideRect(ex2,ey2,rx,ry,rw,rh);
                if (in1 && in2) remove_element_by_2nodes("D", ex1, ey1, ex2, ey2);
                else keep.push_back(d);
            }
            diodes.swap(keep);
        }
        { // V: از x تا x+LEN_V()
            std::vector<Voltage> keep;
            for (auto &v : voltages) {
                int ex1=v.x, ey1=v.y, ex2=v.x + LEN_V(), ey2=v.y;   // ← به‌جای 70
                bool in1=insideRect(ex1,ey1,rx,ry,rw,rh), in2=insideRect(ex2,ey2,rx,ry,rw,rh);
                if (in1 && in2) remove_element_by_2nodes("V", ex1, ey1, ex2, ey2);
                else keep.push_back(v);
            }
            voltages.swap(keep);
        }
        { // I: از x تا x+LEN_I()
            std::vector<CurrentSrc> keep;
            for (auto &i : currents) {
                int ex1=i.x, ey1=i.y, ex2=i.x + LEN_I(), ey2=i.y;   // ← به‌جای 70
                bool in1=insideRect(ex1,ey1,rx,ry,rw,rh), in2=insideRect(ex2,ey2,rx,ry,rw,rh);
                if (in1 && in2) remove_element_by_2nodes("I", ex1, ey1, ex2, ey2);
                else keep.push_back(i);
            }
            currents.swap(keep);
        }
        { // L: از x تا x+LEN_L()
            std::vector<Inductor> keep;
            for (auto &ind : inductors) {
                int ex1=ind.x, ey1=ind.y, ex2=ind.x + LEN_L(), ey2=ind.y;   // ← به‌جای 108
                bool in1=insideRect(ex1,ey1,rx,ry,rw,rh), in2=insideRect(ex2,ey2,rx,ry,rw,rh);
                if (in1 && in2) remove_element_by_2nodes("L", ex1, ey1, ex2, ey2);
                else keep.push_back(ind);
            }
            inductors.swap(keep);
        }

        { // GND: آیکون 28x28 باید کامل داخل مستطیل باشد
            std::vector<Ground> keep;
            for (auto &g : grounds) {
                bool fullyInside = (g.x >= rx && g.x+28 <= rx+rw && g.y >= ry && g.y+28 <= ry+rh);
                if (!fullyInside) keep.push_back(g);
                else {
                    int gid = getExistingNodeIdAt(g.x+14, g.y);
                    if (gid >= 0) {
                        for (size_t j=0; j<ground_node_ids.size(); ) {
                            if (ground_node_ids[j]==gid) ground_node_ids.erase(ground_node_ids.begin()+j);
                            else ++j;
                        }
                    }
                }
            }
            grounds.swap(keep);

        }

        // B label: نقطه 5x5 باید "کامل" داخل مستطیل باشد
        {
            std::vector<LabelB> keep;
            for (auto &b : labelsB) {
                int left   = b.x - 2;
                int top    = b.y - 2;
                int right  = left + 5;
                int bottom = top  + 5;
                bool fullyInside = (left >= rx && right <= rx+rw && top >= ry && bottom <= ry+rh);
                if (!fullyInside) keep.push_back(b);
                // اگر fullyInside بود، عمداً نگه نمی‌داریم ⇒ حذف می‌شود
            }
            labelsB.swap(keep);


        }
        cleanupDanglingNodes();

        // 4) سینک کردن elements با wires باقی‌مانده (حذف R با مقدار سیم)
        {
            set<std::pair<int,int>> wirePairs;
            for (auto &w : wires) {
                int id1 = getExistingNodeIdAt(w.x1, w.y1);
                int id2 = getExistingNodeIdAt(w.x2, w.y2);
                if (id1>=0 && id2>=0) {
                    if (id1>id2) std::swap(id1,id2);
                    wirePairs.insert({id1,id2});
                }
            }
            for (size_t i=0; i<elements.size(); ) {
                if (elements[i].type=="R" && std::fabs(elements[i].value - WIRE_RESISTANCE) < 1e-9) {
                    int a = std::min(elements[i].node1, elements[i].node2);
                    int b = std::max(elements[i].node1, elements[i].node2);
                    if (!wirePairs.count({a,b})) { elements.erase(elements.begin()+i); continue; }
                }
                ++i;
            }
        }

        // 5) پاکسازی گره‌های بی‌مصرف و بازسازی map/parent (مثل قبلی)
        {
            std::vector<char> used(nodes.size(), 0);
            for (auto &el : elements) { used[el.node1]=1; used[el.node2]=1; }
            for (int gid : ground_node_ids) if (gid>=0 && gid<(int)used.size()) used[gid]=1;

            for (auto &b : labelsB)
                if (b.nodeId >= 0 && b.nodeId < (int)used.size())
                    used[b.nodeId] = 1;

            std::vector<int> mapOldToNew(nodes.size(), -1);
            std::vector<Node> newNodes;
            for (int i=0; i<(int)nodes.size(); ++i) if (used[i]) {
                mapOldToNew[i] = (int)newNodes.size();
                Node nn = nodes[i]; nn.id = (int)newNodes.size(); newNodes.push_back(nn);
            }
            for (auto &el : elements) { el.node1 = mapOldToNew[el.node1]; el.node2 = mapOldToNew[el.node2]; }
            for (size_t i=0; i<ground_node_ids.size(); ) {
                int old = ground_node_ids[i];
                int neu = (old>=0 && old<(int)mapOldToNew.size()) ? mapOldToNew[old] : -1;
                if (neu==-1) ground_node_ids.erase(ground_node_ids.begin()+i); else { ground_node_ids[i]=neu; ++i; }
            }

            for (auto &b : labelsB) {
                int old = b.nodeId;
                int neu = (old >= 0 && old < (int)mapOldToNew.size()) ? mapOldToNew[old] : -1;
                if (neu == -1) {
                    // اگر به هر دلیلی گره حذف شده بود، Label را هم بندازیم
                    b.nodeId = -1;
                } else {
                    b.nodeId = neu;
                }
            }
            // (اختیاری) پاک کردن Labelهایی که nodeId شان -1 شد:
            labelsB.erase(std::remove_if(labelsB.begin(), labelsB.end(),
                [](const LabelB& b){ return b.nodeId < 0; }), labelsB.end());

            // پس از حذف برچسب‌ها:
            std::set<int> alive;
            for (auto &b : labelsB) if (b.nodeId >= 0) alive.insert(b.nodeId);

            for (auto &g : labelB_groups) {
                g.erase(std::remove_if(g.begin(), g.end(),
                        [&](int id){ return !alive.count(id); }), g.end());
            }
            labelB_groups.erase(std::remove_if(labelB_groups.begin(), labelB_groups.end(),
                    [](const std::vector<int>& g){ return g.size() < 2; }),
                    labelB_groups.end());

            nodes.swap(newNodes);
            parent.assign(nodes.size(), 0);
            for (int i=0; i<(int)nodes.size(); ++i) parent[i]=i;
        }
    }
}




        }

        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderClear(renderer);
        draw_editor_grid(renderer, SCREEN_WIDTH, SCREEN_HEIGHT, GRID_SIZE);
        rect(renderer, 0, 0, SCREEN_WIDTH, 73, 200, 200, 200, 255, 1);
         texture(renderer, 700, 20, "resistor1.png", 26, 28);
        texture(renderer, 735, 20, "capacitor1.png", 26, 28);
        texture(renderer, 770, 20, "diod.png", 26, 28);
        texture(renderer, 805, 20, "dc-voltage-source.png", 26, 28);
        texture(renderer, 840, 20, "inductor1.png", 26, 28);
        texture(renderer, 875, 20, "gnd.png", 28, 28);
        texture(renderer, 910, 20, "Current_Source.png", 26, 28); // آیکون منبع جریان – بعد از GND

        stringRGBA(renderer, 952, 26, (char*)"B", 0, 0, 0, 255); // حرف B
        SDL_Rect bDotIcon = { 964, 39, 5, 5 };                   // نقطه 5x5 کنار B
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderFillRect(renderer, &bDotIcon);


        texture(renderer, 20, 20, "new1.png", 26, 28);
        texture(renderer, 55, 20, "save.png", 26, 28);
        texture(renderer, 90, 20, "run.png", 26, 28);
        texture(renderer, 665, 20, "cabel.png", 26, 28);
        texture(renderer, 630, 20, "delete.png", 26, 28); // آیکون حذف انتخابی

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        for (auto& r : resistors) {
            draw_resistor(renderer, r.x, r.y, r.angle);
        }
        for (auto& c : capacitors) draw_capacity(renderer, c.x, c.y);
        for (auto& d : diodes) draw_diode(renderer, d.x, d.y);
        for (auto& v : voltages) draw_voltage(renderer, v.x, v.y);
        for (auto& i : currents) draw_current_source(renderer, i.x, i.y);
        for (auto& i : inductors) draw_inductor(renderer, i.x, i.y);
        for (auto& g : grounds) draw_ground(renderer, g.x, g.y);
        for (auto &b : labelsB) {
            SDL_Rect s = { b.x - 2, b.y - 2, 5, 5 };        // نقطه 5x5
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
            SDL_RenderFillRect(renderer, &s);
            stringRGBA(renderer, b.x + 8, b.y - 8, (char*)"B", 0, 0, 0, 255); // متن B
        }
        for (auto& w : wires) {
            SDL_RenderDrawLine(renderer, w.x1, w.y1, w.x2, w.y2);
            draw_connector(renderer, w.x1, w.y1);
            draw_connector(renderer, w.x2, w.y2);
        }
        if (showNodeIds) {
            for (const auto& n : nodes) {
                char buf[16];
                sprintf(buf, "%d", n.id);
                // کمی راست-بالا نسبت به خود گره تا با مربع قرمز تداخل نداشته باشه
                stringRGBA(renderer, n.x + 6, n.y - 10, (char*)buf, 0, 0, 0, 255);
            }
        }

        int mx, my; SDL_GetMouseState(&mx, &my);
        mx = snap(mx); my = snap(my);

        if (cursorMode == PLACING_RESISTOR) {
            // جدید (Anchor = مرکز ماوس → چپ = مرکز - LEN_R()/2):
            int left = snap(mx - LEN_R()/2);
            draw_resistor(renderer, left, my, current_rotation);
        }   else if (cursorMode == PLACING_CAPACITOR) {
draw_capacity(renderer,        mx - LEN_C()/2, my);
        } else if (cursorMode == PLACING_DIODE) {
            draw_diode(renderer,           mx - LEN_D()/2, my);

        }  else if (cursorMode == PLACING_LABEL_B) {
            SDL_Rect s = { mx - 2, my - 2, 5, 5 };
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
            SDL_RenderFillRect(renderer, &s);
            stringRGBA(renderer, mx + 8, my - 8, (char*)"B", 0, 0, 0, 255);
        }else if (cursorMode == PLACING_VOLTAGE) {
            draw_voltage(renderer,         mx - LEN_V()/2, my);
        } else if (cursorMode == PLACING_CURRENT_SOURCE) {
            draw_current_source(renderer,  mx - LEN_I()/2, my);
        }else if (cursorMode == PLACING_INDUCTOR) {
            draw_inductor(renderer,        mx - LEN_L()/2, my);
        }
        else if (cursorMode == PLACING_GROUND) {
            draw_ground(renderer, mx - 14, my);
        } else if (cursorMode == DRAWING_WIRE && wire_in_progress && preview_has) {
            SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255); // اگر شفاف می‌خواهی: BlendMode را فعال کن
            SDL_RenderDrawLine(renderer, wire_start_x, wire_start_y, preview_x2, preview_y2);
            // (اختیاری) کانکتورهای دو سر پیش‌نمایش:
             draw_connector(renderer, wire_start_x, wire_start_y);
             draw_connector(renderer, preview_x2, preview_y2);
        }

        if (cursorMode == DELETE_SELECT && del_select_active) {
            int rx, ry, rw, rh;
            normRect(sel_x0, sel_y0, sel_x1, sel_y1, rx, ry, rw, rh);
            SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
            SDL_Rect r{rx, ry, rw, rh};
            SDL_RenderDrawRect(renderer, &r);
        }


        SDL_RenderPresent(renderer);
    }

    SDL_FreeCursor(cursor);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    TTF_Quit();
    SDL_Quit();

    return 0;
}

void rect(SDL_Renderer* Renderer, int x, int y, int w, int h, int R, int G, int B, int A, int fill_boolian) {
    SDL_Rect rectangle = { x, y, w, h };
    SDL_SetRenderDrawColor(Renderer, R, G, B, A);
    if (fill_boolian)
        SDL_RenderFillRect(Renderer, &rectangle);
    SDL_RenderDrawRect(Renderer, &rectangle);
}

void texture(SDL_Renderer* renderer, int x, int y, string address, int w, int h) {
    int n = address.length();
    char chhar_arr[n + 1];
    strcpy(chhar_arr, address.c_str());
    SDL_Texture* myTexture = IMG_LoadTexture(renderer, chhar_arr);
    SDL_Rect texr1 = { x, y, w, h };
    SDL_RenderCopy(renderer, myTexture, nullptr, &texr1);
    SDL_DestroyTexture(myTexture);
}

void rotate_point(int cx, int cy, int& x, int& y, int angle) {
    double rad = angle * M_PI / 180.0;
    int tx = x - cx;
    int ty = y - cy;
    int rx = static_cast<int>(tx * cos(rad) - ty * sin(rad));
    int ry = static_cast<int>(tx * sin(rad) + ty * cos(rad));
    x = cx + rx;
    y = cy + ry;
}



void draw_ground(SDL_Renderer* renderer, int x, int y) {
    // Draw 28x28 ground icon. User said the file is already placed next to others.
    texture(renderer, x, y, "gnd.png", 28, 28);
    // Single 5x5 connector at the top-center of the icon
    draw_connector(renderer, x + 14, y);
}

void draw_resistor(SDL_Renderer* renderer, int x, int y, int angle) {
    const int L = LEN_R();
    const double step = double(L) / 7.0;  // 8 segment points (0..7)
    int amp = GRID_SIZE / 2; if (amp < 6) amp = 6;

    int px[8], py[8];
    for (int i = 0; i < 8; ++i) {
        px[i] = x + int(std::round(i * step));
        py[i] = y;
    }
    // زیگزاگ
    py[2] = y - amp; py[3] = y + amp; py[4] = y - amp; py[5] = y + amp;

    // چرخش حول مرکز
    int cx = x + L/2, cy = y;
    for (int i = 0; i < 8; ++i) rotate_point(cx, cy, px[i], py[i], angle);

    // رسم
    for (int i = 0; i < 7; ++i) SDL_RenderDrawLine(renderer, px[i], py[i], px[i+1], py[i+1]);
    draw_connector(renderer, px[0], py[0]);
    draw_connector(renderer, px[7], py[7]);
}



void draw_capacity(SDL_Renderer* renderer, int x, int y, int angle) {
    const int L = LEN_C();
    const int gap = GRID_SIZE;           // فاصله بین دو صفحه
    const int lead = (L - gap) / 2;      // طول سیم‌های دو طرف
    const int plateH = GRID_SIZE;        // ارتفاع صفحات

    int x1 = x + lead;        // صفحهٔ اول
    int x2 = x1 + gap;        // صفحهٔ دوم

    // (در این نسخه چرخش اعمال نمی‌کنیم؛ اگر زاویه لازم داشتی مشابه مقاومت حول مرکز بچرخان)
    SDL_RenderDrawLine(renderer, x, y, x1, y);                  // سیم چپ
    SDL_RenderDrawLine(renderer, x2, y, x + L, y);              // سیم راست
    SDL_RenderDrawLine(renderer, x1, y - plateH/2, x1, y + plateH/2);
    SDL_RenderDrawLine(renderer, x2, y - plateH/2, x2, y + plateH/2);

    draw_connector(renderer, x, y);
    draw_connector(renderer, x + L, y);
}


void draw_diode(SDL_Renderer* renderer, int x, int y, int angle) {
    const int L = LEN_D();
    const int u = L / 3;                 // lead / body / lead
    const int h = GRID_SIZE / 2; if (h < 6) {} // ارتفاع

    // سیم‌های چپ و راست
    SDL_RenderDrawLine(renderer, x, y, x + u, y);
    SDL_RenderDrawLine(renderer, x + 2*u, y, x + L, y);

    // بدنهٔ دیود (مثلث + خط کاتد در x+2u)
    SDL_RenderDrawLine(renderer, x + u, y - h, x + u, y + h);          // ضلع عمودی مثلث
    SDL_RenderDrawLine(renderer, x + u, y - h, x + 2*u, y);            // ضلع بالا
    SDL_RenderDrawLine(renderer, x + u, y + h, x + 2*u, y);            // ضلع پایین
    SDL_RenderDrawLine(renderer, x + 2*u, y - h, x + 2*u, y + h);      // کاتد

    draw_connector(renderer, x, y);
    draw_connector(renderer, x + L, y);
}


void draw_voltage(SDL_Renderer* renderer, int x, int y, int angle) {
    const int L = LEN_V();
    const int r = L / 4;  // شعاع دایره
    // سیم‌ها
    SDL_RenderDrawLine(renderer, x, y, x + r, y);
    SDL_RenderDrawLine(renderer, x + 3*r, y, x + L, y);
    // دایره
    circleRGBA(renderer, x + L/2, y, r, 0, 0, 0, 255);
    // علامت + و - داخل دایره (ساده)
    SDL_RenderDrawLine(renderer, x + L/2 - 8, y, x + L/2 - 18, y);
    SDL_RenderDrawLine(renderer, x + L/2 + 8, y - 5, x + L/2 + 8, y + 5);

    draw_connector(renderer, x, y);
    draw_connector(renderer, x + L, y);
}

void draw_current_source(SDL_Renderer* renderer, int x, int y, int angle) {
    const int L = LEN_I();
    const int r = L / 4;
    SDL_RenderDrawLine(renderer, x, y, x + r, y);
    SDL_RenderDrawLine(renderer, x + 3*r, y, x + L, y);
    circleRGBA(renderer, x + L/2, y, r, 0, 0, 0, 255);
    // فلش داخل دایره
    int ax1 = x + L/2 - 8, ax2 = x + L/2 + 8;
    SDL_RenderDrawLine(renderer, ax1, y, ax2, y);
    SDL_RenderDrawLine(renderer, ax2, y, ax2 - 5, y - 5);
    SDL_RenderDrawLine(renderer, ax2, y, ax2 - 5, y + 5);

    draw_connector(renderer, x, y);
    draw_connector(renderer, x + L, y);
}
void draw_inductor(SDL_Renderer* renderer, int x, int y, int angle) {
    const int L = LEN_L();
    const int radius  = GRID_SIZE / 2;     // 10 (با GRID_SIZE=20)
    const int spacing = GRID_SIZE / 2;     // 10
    const int count   = 3;                 // سه نیم‌دایره
    const int coilsW  = count * 2 * radius;
    const int leftLead = spacing;
    const int rightLead = spacing;

    // سرها
    draw_connector(renderer, x, y);
    draw_connector(renderer, x + L, y);

    // سیم‌های صاف
    SDL_RenderDrawLine(renderer, x, y, x + leftLead, y);
    SDL_RenderDrawLine(renderer, x + leftLead + coilsW, y, x + leftLead + coilsW + rightLead, y);

    // حلقه‌ها
    int startx = x + leftLead;
    for (int i = 0; i < count; ++i) {
        for (int a = 0; a <= 180; a += 10) {
            float rad = a * float(M_PI) / 180.0f;
            int px = startx + i * (2*radius) + radius + std::cos(rad) * radius;
            int py = y + std::sin(rad) * radius;
            SDL_RenderDrawPoint(renderer, px, py);
        }
    }
}

// قبلی را با این عوض کن
void exportToFile(const std::string& path){
    // Ensure merged nodes share same parent (همان کد union-find فعلی)
    for (int i = 0; i < (int)nodes.size(); i++) parent[i] = findParent(i);

    std::map<int,int> nodeMap;
    std::map<int,int> parentToNum;
    int nextNodeNum = 1;
    for (int i = 0; i < (int)nodes.size(); i++) {
        int p = parent[i];
        if (!parentToNum.count(p)) parentToNum[p] = nextNodeNum++;
        nodeMap[i] = parentToNum[p];
    }

    int rCount = 1, cCount = 1, lCount = 1, vCount = 1, iCount = 1, dCount = 1;
    int rbCount=1; // شمارنده‌ی جدا برای سیم نامرئی

    std::ofstream out(path);

    // --- GND export: if any ground placed, take the last one as the reference.
    if (!ground_node_ids.empty()) {
        int lastGroundNode = ground_node_ids.back();
        int mapped = nodeMap[lastGroundNode];
        out << "NODE GND " << mapped << "\n";
    }

    if (!out.is_open()) {
        std::cerr << "Error: cannot open " << path << " for writing\n";
        return;
    }

    for (size_t idx = 0; idx < elements.size(); ++idx) {
        auto &el = elements[idx];
        int n1 = nodeMap[el.node1];
        int n2 = nodeMap[el.node2];

        if (el.type == "R") {
            out << "R R" << rCount++ << " " << n1 << " " << n2 << " " << el.value << "\n";
        } else if (el.type == "C") {
            out << "C C" << cCount++ << " " << n1 << " " << n2 << " " << el.value << "\n";
        } else if (el.type == "L") {
            out << "L L" << lCount++ << " " << n1 << " " << n2 << " " << el.value << "\n";
        } else if (el.type == "V") {
            auto it = srcWFByIdx.find((int)idx);
            if (it != srcWFByIdx.end()) {
                const auto &wf = it->second;
                if (wf.isSin) {
                    out << "V V" << vCount++ << " " << n1 << " " << n2
                        << " " << "SIN(" << wf.off << " " << wf.amp << " " << wf.freq << ")\n";
                } else if (wf.isPulse) {
                    out << "V V" << vCount++ << " " << n1 << " " << n2
                        << " " << "PULSE " << wf.pulseType << " (" << wf.off << " " << wf.amp << " " << wf.freq << ")\n";
                } else {
                    out << "V V" << vCount++ << " " << n1 << " " << n2 << " " << el.value << "\n";
                }
            } else {
                out << "V V" << vCount++ << " " << n1 << " " << n2 << " " << el.value << "\n";
            }

        } else if (el.type == "I") {
            auto it = srcWFByIdx.find((int)idx);
            if (it != srcWFByIdx.end()) {
                const auto &wf = it->second;
                if (wf.isSin) {
                    out << "I I" << iCount++ << " " << n1 << " " << n2
                        << " " << "SIN(" << wf.off << " " << wf.amp << " " << wf.freq << ")\n";
                } else if (wf.isPulse) {
                    out << "I I" << iCount++ << " " << n1 << " " << n2
                        << " " << "PULSE " << wf.pulseType << " (" << wf.off << " " << wf.amp << " " << wf.freq << ")\n";
                } else {
                    out << "I I" << iCount++ << " " << n1 << " " << n2 << " " << el.value << "\n";
                }
            } else {
                out << "I I" << iCount++ << " " << n1 << " " << n2 << " " << el.value << "\n";
            }

        } else if (el.type == "D") {
            out << "D D" << dCount++ << " " << n1 << " " << n2 << " D\n";
        }
        for (const auto& grp : labelB_groups) {             // ← توضیح بخش 2 را پایین ببینید
            for (size_t k=1; k<grp.size(); ++k) {
                int n1 = nodeMap[grp[k-1]];
                int n2 = nodeMap[grp[k]];
                if (n1!=n2)
                    out << "R RB" << rbCount++ << " " << n1 << " " << n2 << " " << WIRE_RESISTANCE << "\n";
            }
        }
    }

    out.close();
    std::cout << "Exported circuit to " << path << "\n";
}



// --- NEW: پاک‌سازی کامل صفحه و ریست وضعیت ---
static void clearAll() {
        resistors.clear();
        capacitors.clear();
        diodes.clear();
        voltages.clear();
        currents.clear();          // اگر نیست، با پچ قبلی اضافه شده است
        inductors.clear();
        grounds.clear();
        ground_node_ids.clear();
        wires.clear();
        nodes.clear();
        elements.clear();
        parent.clear();
    labelsB.clear();
    labelB_session_count = 0;
        wire_in_progress = false;
        del_select_active = false;
    }

// --- NEW: گرفتن مسیر ذخیره با پنجرهٔ Save ویندوز ---
static std::string ask_user_for_save_path() {
        OPENFILENAMEA ofn{};
        char fileBuf[MAX_PATH]   = {0};
        char initDir[MAX_PATH]   = {0};
        char originalDir[MAX_PATH] = {0};
        GetCurrentDirectoryA(MAX_PATH, initDir);
        GetCurrentDirectoryA(MAX_PATH, originalDir); // مسیر فعلی را نگه دار

        ofn.lStructSize     = sizeof(ofn);
        ofn.hwndOwner       = nullptr;
        ofn.lpstrFile       = fileBuf;
        ofn.nMaxFile        = MAX_PATH;
        ofn.lpstrFilter     = "Text files (*.txt)\0*.txt\0All files (*.*)\0*.*\0";
        ofn.nFilterIndex    = 1;
        ofn.lpstrInitialDir = initDir;
        ofn.Flags           = OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST | OFN_NOCHANGEDIR; // مسیر جاری را تغییر نده
        ofn.lpstrDefExt     = "txt";

        std::string path;
        if (GetSaveFileNameA(&ofn)) path = std::string(ofn.lpstrFile);

        // هر نتیجه‌ای که بود، مسیر جاری برنامه را به حالت قبل برگردان:
        SetCurrentDirectoryA(originalDir);
        return path;
    }


static void draw_editor_grid(SDL_Renderer* r, int w, int h, int spacing) {
    // خطوط کم‌رنگ
    SDL_SetRenderDrawColor(r, 230,230,230,255);
    for (int x = 0; x <= w; x += spacing) SDL_RenderDrawLine(r, x, 0, x, h);
    for (int y = 0; y <= h; y += spacing) SDL_RenderDrawLine(r, 0, y, w, y);
    // هر 5 خط را پررنگ‌تر کن
    SDL_SetRenderDrawColor(r, 210,210,210,255);
    for (int x = 0, k=0; x <= w; x += spacing,++k) if (k%5==0) SDL_RenderDrawLine(r, x, 0, x, h);
    for (int y = 0, k=0; y <= h; y += spacing,++k) if (k%5==0) SDL_RenderDrawLine(r, 0, y, w, y);
}


static bool snapToNearestNode(int x, int y, int& sx, int& sy) {
    int best = -1, bestd2 = MAGNET_RADIUS*MAGNET_RADIUS;
    for (int i=0;i<(int)nodes.size();++i) {
        int dx = nodes[i].x - x, dy = nodes[i].y - y;
        int d2 = dx*dx + dy*dy;
        if (d2 <= bestd2) { bestd2 = d2; best = i; }
    }
    if (best >= 0) { sx = nodes[best].x; sy = nodes[best].y; return true; }
    return false;
}


static int pickNodeAt(int mx, int my) {
    for (int i=0;i<(int)nodes.size();++i) {
        int dx = nodes[i].x - mx, dy = nodes[i].y - my;
        if (dx*dx + dy*dy <= NODE_HIT_RADIUS*NODE_HIT_RADIUS) return i;
    }
    return -1;
}
// نام‌نمایشی: R1, C2, ...
static std::string elName(size_t idx) {
    if (idx >= elements.size()) return "?";
    const Element& el = elements[idx];
    int count = 0;
    for (size_t i = 0; i <= idx; ++i)
        if (elements[i].type == el.type) ++count;
    return el.type + std::to_string(count);
}


static void cleanupDanglingNodes() {
    // 1) علامت‌گذاری نودهای استفاده‌شده
    std::vector<char> used(nodes.size(), 0);
    for (auto &el : elements) {
        if (el.node1 >= 0 && el.node1 < (int)used.size()) used[el.node1] = 1;
        if (el.node2 >= 0 && el.node2 < (int)used.size()) used[el.node2] = 1;
    }
    for (int g : ground_node_ids)
        if (g >= 0 && g < (int)used.size()) used[g] = 1;
    for (auto &b : labelsB)
        if (b.nodeId >= 0 && b.nodeId < (int)used.size()) used[b.nodeId] = 1;

    // 2) ساخت مپ قدیم→جدید و کپی نودهای زنده
    std::vector<int> mapOldToNew(nodes.size(), -1);
    std::vector<Node> newNodes; newNodes.reserve(nodes.size());
    for (int i = 0; i < (int)nodes.size(); ++i) if (used[i]) {
        mapOldToNew[i] = (int)newNodes.size();
        Node nn = nodes[i]; nn.id = (int)newNodes.size();
        newNodes.push_back(nn);
    }

    // 3) ری‌مپ المنت‌ها
    for (auto &el : elements) {
        el.node1 = mapOldToNew[el.node1];
        el.node2 = mapOldToNew[el.node2];
    }

    // 4) ری‌مپ گراندها (و حذف گراندهای بی‌اعتبار)
    for (size_t i = 0; i < ground_node_ids.size(); ) {
        int old = ground_node_ids[i];
        int neu = (old >= 0 && old < (int)mapOldToNew.size()) ? mapOldToNew[old] : -1;
        if (neu == -1) ground_node_ids.erase(ground_node_ids.begin()+i);
        else { ground_node_ids[i] = neu; ++i; }
    }

    // 5) ری‌مپ/حذف LabelBهای آویزان و تمیز کردن گروه‌ها
    for (auto &b : labelsB) {
        int old = b.nodeId;
        int neu = (old >= 0 && old < (int)mapOldToNew.size()) ? mapOldToNew[old] : -1;
        b.nodeId = neu;
    }
    labelsB.erase(std::remove_if(labelsB.begin(), labelsB.end(),
                  [](const LabelB& b){ return b.nodeId < 0; }), labelsB.end());

    std::set<int> alive;
    for (auto &b : labelsB) if (b.nodeId >= 0) alive.insert(b.nodeId);
    for (auto &g : labelB_groups) {
        g.erase(std::remove_if(g.begin(), g.end(),
               [&](int id){ return !alive.count(id); }), g.end());
    }
    labelB_groups.erase(std::remove_if(labelB_groups.begin(), labelB_groups.end(),
               [](const std::vector<int>& g){ return g.size() < 2; }),
               labelB_groups.end());

    // 6) جاگذاری نودهای جدید و ریست union-find
    nodes.swap(newNodes);
    parent.assign(nodes.size(), 0);
    for (int i = 0; i < (int)nodes.size(); ++i) parent[i] = i;
}
