# راهنمای آپلود فایل‌های پروژه در Repository دانشگاه

این پروژه قرار است داخل **همان مخزن ایجادشده توسط استاد** نگهداری شود؛ لازم نیست Repository جدید بسازید.

## اولین بار

1. آدرس دقیق Repository دانشگاه را بردارید؛ دسترسی و سیاست‌های دانشگاه (SSO یا SSH) را چک کنید.
2. در PowerShell یا Git Bash دستور زیر را با URL واقعی اجرا کنید:

```bash
git clone <UNIVERSITY_REPOSITORY_URL>
cd <FOLDER_NAME>
git status
```

3. محتوای ZIP اولیه را **داخل ریشهٔ پوشهٔ مخزن Clone شده** کپی کنید. اگر آنجا README یا .gitignore وجود دارد، آنها را **ادغام** کنید، نه اینکه رویشان بدون بررسی بنویسید.
4. مطابق قواعد استاد یک Branch بسازید و فایل‌های متنی و کدها را ثبت کنید:

```bash
git switch -c feature/step-01-project-scope
git add docs/ config/ data/metadata/ notebooks/ src/ tests/ results/ thesis/ environment.yml
git status
git diff --cached --stat
git commit -m "docs: initialize thesis scope and repository layout"
git push -u origin feature/step-01-project-scope
```

5. در صورت استفادهٔ گروه از Pull Request، تغییرات را برای تأیید استاد ارسال کنید. README و .gitignore را پس از بررسی نسخهٔ موجود اضافه کنید.

## فایل‌ها کجا می‌روند؟
- کد Python: `src/`
- Jupyter Notebook: `notebooks/`
- روش، فرضیات، مقاله‌ها و تصمیم‌ها: `docs/`
- پیکربندی: `config/`
- فهرست منبع و مجوز داده‌ها: `data/metadata/sources.csv`
- نقشه‌ها/نمودارهای **منتخب و کوچک**: `results/figures/`
- جدول‌های منتخب: `results/tables/`
- فایل‌های خود پایان‌نامه: `thesis/`

## چه چیزهایی را آپلود نکنیم؟
داده‌های خام حجیم `.tif`, `.gpkg`, `.shp`، فایل‌های محدودشدهٔ دانشگاه یا اشخاص، رمزها و توکن‌ها، خروجی‌های موقت و کش‌ها. برای اشتراک داده‌ها از فضای مورد تأیید دانشگاه استفاده کنید و محل دسترسی و مجوز را در manifest ثبت کنید.

## بعد از هر استپ

```bash
git status
git add docs/ src/ notebooks/ config/ data/metadata/ results/figures/ results/tables/
git diff --cached --stat
git commit -m "docs: update study methods"
git push
```

**توجه:** اسم Branch، درخواست Pull Request و محل نگهداری داده باید با قواعد Repository استاد هماهنگ شود.
