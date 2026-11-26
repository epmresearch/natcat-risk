# راهنمای فارسی اجرای Step 4.1 — فهرست منابع داده و کنترل محدودهٔ مطالعه

## هدف

در این مرحله **نه مدل خطر ساخته می‌شود و نه DEM حجیم دانلود می‌کنیم**. منابع قابل استناد شناسایی و ثبت می‌شوند و اسکریپت جدید بررسی می‌کند که فایل محدودهٔ جست‌وجو و فهرست منابع قابل استفاده و سازگار هستند. هنوز هیچ‌یک از منابع بیرونی برای پوشش تک‌تک نقاط Highway 99 بررسی مکانی نشده‌اند.

## فایل‌های پکیج

- `data/catalog/step4_1_source_registry.csv` — نه منبع پیشنهادی، همراه با اولویت، پیوند، کاربرد، وضعیت مجوز و محدودیت‌ها.
- `docs/step4_1_data_source_assessment.md` — یادداشت فنی دانشگاهی انگلیسی با منابع و معیارهای انتخاب.
- `docs/step4_1_execution_FA.md` — همین راهنما.
- `src/data_sources/01_audit_step41_sources.py` — کنترل آفلاین لیست منابع و محدودهٔ ۵ کیلومتری.
- `tests/test_step4_1_data_sources.py` — ۱۱ تست مستقل کد.

## اجرای محلی

در **Miniforge Prompt**:

```bat
mamba activate sea_to_sky_risk
cd /d "F:\Projects\University of Calgary\Thesis\natcat-risk"
python -m unittest discover -s tests -p "test_step4_1_data_sources.py"
python src/data_sources/01_audit_step41_sources.py
```

خروجی دوم در مسیر زیر ایجاد می‌شود:

`results/diagnostics/step4_1_source_discovery_audit.json`

## انتظار از خروجی

- شمار منابع: `9`.
- مساحت محدودهٔ ۵ کیلومتری: حدود `1043.616776 km²`.
- تطابق مساحت با گزارش Step 3.4: `true`.
- وضعیت: `SOURCE_REGISTRY_AND_AOI_QA_PASSED_COVERAGE_NOT_VERIFIED`.
- شمار منابعی که پوشش محلی‌شان *مستقل بررسی شده*: `0` (در این مرحله کاملاً طبیعی است).

## بعد از اجرا چه بفرستی؟

۱. خروجی کامل Terminal؛ ۲. فایل `results/diagnostics/step4_1_source_discovery_audit.json`.

پس از بررسی، **داخل QGIS یا با سرویس STAC** پوشش دادهٔ ارتفاعی MRDEM-30 و HRDEM رو واقعاً کنترل می‌کنیم. راهنمای رسمی QGIS سازمان منابع طبیعی کانادا: https://nrcan.github.io/CanElevation/stac-dem-mosaics/accessing-stac-with-qgis/ .

تا زمانی که لینک‌ها، پوشش واقعی، سیستم ارتفاع و مقیاس بررسی نشوند، دانلود انبوه نکن. هنوز Git Commit هم نزن؛ ابتدا خروجی کنترل رو بررسی می‌کنیم.

## مرز روش علمی

- ۵ کیلومتر بافر یعنی **محدودهٔ جست‌وجوی داده**، نه محدودهٔ خطر.
- DEM ۳۰ متری به معنی مشاهدهٔ دقیق شیب برش جاده نیست.
- پایگاه v14 زمین‌لغزش شامل عوارض و رخدادها در کل کانادا است؛ بدون فیلتر و کنترل زمانی نمی‌توان آمار وقوع در Highway 99 از آن نتیجه گرفت.
- مرکزخط Highway 99 هنوز OSM-based و *provisional* است.
